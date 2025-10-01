module Experiments

export Experiment, ExperimentCfg, CheckpointCfg, run_experiment, load_experiment
export DelayedCfg, ManyColsCfg
export load_data

using Serialization, JSON3, Setfield, DataStructures, Dates

using Mavi.Systems
using Mavi.MaviSerder
using Mavi.Integration: get_step_function

using Mavi.Utils.Progress

# Module Configs
COL_DIRNAME = "data"
FINAL_STATE_NAME = "final_state"
CHECKPOINT_DIRNAME = "checkpoint"
CHECKPOINT_INFO_NAME = "checkpoint_info.bin"
EXP_COL_CONFIGS_NAME = "experiment_configs.json"

abstract type ColCfg end
abstract type ColState end
abstract type Collector end

function get_collector(col_cfg::ColCfg, exp_cfg, system, state=nothing) end
function collect(col::Collector, system) end
function save_data(col::Collector, path) end
function load_data(::Type{ColCfg}, path) end

# = 
# Experiment
# = 

@kwdef struct CheckpointCfg
    delta_time::Float64
    save_data::Bool=false
end

mutable struct Checkpoint
    last_used::String
    last_time::Float64
    is_valid::Dict{String, Bool}
end
Checkpoint(last_time) = Checkpoint("2",  last_time, Dict("1"=>false, "2"=>false))

@kwdef struct ExperimentCfg 
    tf::Float64
    root::String
    save_final_state::Bool=true
    checkpoint_cfg::Union{Nothing, CheckpointCfg}=nothing
end

struct Experiment{S<:System, F<:Union{Nothing, Function}, C<:Union{Nothing, Checkpoint}} 
    cfg::ExperimentCfg
    col::Collector
    system::S
    custom_step::F
    checkpoint::C
end
function Experiment(
    cfg::ExperimentCfg, col_cfg::ColCfg, system::System, 
    custom_step=nothing)
    
    if isnothing(custom_step)
        custom_step = get_step_function(system.type, system)
    end

    checkpoint = nothing
    if !isnothing(cfg.checkpoint_cfg)
        checkpoint = Checkpoint(system.time_info.time)
    end

    Experiment(cfg, get_collector(col_cfg, cfg, system), system, custom_step, checkpoint)
end

function run_experiment(experiment::Experiment, stop_func=nothing)
    system = experiment.system
    col = experiment.col
    cfg = experiment.cfg

    if isnothing(stop_func)
        stop_func = (system) -> false
    end

    save_system_configs(system, cfg.root)

    open(joinpath(cfg.root, EXP_COL_CONFIGS_NAME), "w") do io
        JSON3.pretty(io, 
            Dict(
                "experiment" => get_obj_save_data_json(cfg), 
                "collector" => get_obj_save_data_json(col.cfg), 
            ),
        )
    end

    col_path = mkpath(joinpath(cfg.root, COL_DIRNAME))

    if !isnothing(experiment.checkpoint)
        mkpath(joinpath(cfg.root, CHECKPOINT_DIRNAME, "1"))
        mkpath(joinpath(cfg.root, CHECKPOINT_DIRNAME, "2"))
    end

    experiment_step! = experiment.custom_step

    prog = ProgContinuos(init=system.time_info.time, final=cfg.tf)
    while system.time_info.time < cfg.tf
        experiment_step!(system)
        collect(col, system)
        check_checkpoint(cfg.checkpoint_cfg, experiment)
        show_progress(prog, system.time_info.time)

        if stop_func(system)
            break
        end
    end
    save_data(col, col_path)
    
    if cfg.save_final_state
        save_system(system, joinpath(cfg.root, FINAL_STATE_NAME))
        # save_component_serial(system.state, cfg.root, FINAL_STATE_NAME)
    end

    show_finish(prog)

    open(joinpath(cfg.root, ".done"), "w") do f
        write(f, "Experiment concluded on $(Dates.now())\n")
    end
    
    return
end

# = 
# Collectors Manager
# = 

@kwdef struct ManyColsCfg{C<:NamedTuple} <: ColCfg
    cfgs::C
end

struct ManyColsState{S<:NamedTuple} <: ColState
    states::S
end

struct ManyCols{C, S, ColT<:NamedTuple} <: Collector
    cfg::ManyColsCfg{C}
    state::ManyColsState{S}
    cols::ColT
end

function get_collector(cfg::ManyColsCfg, exp_cfg::ExperimentCfg, system::System, state=nothing)
    cols = []
    cols_states = []
    for (col_name, col_cfg) in pairs(cfg.cfgs)
        col_state = nothing
        if !isnothing(state)
            col_state = state.states[col_name]
        end
        
        col = get_collector(col_cfg, exp_cfg, system, col_state)
        push!(cols, (col_name, col))
        push!(cols_states, (col_name, col.state))
    end

    ManyCols(cfg, ManyColsState(NamedTuple(cols_states)), NamedTuple(cols))
end

function collect(col::ManyCols, system::System)
    for col in values(col.cols)
        collect(col, system)
    end
end

function save_data(col::ManyCols, path)
    for (col_name, col) in pairs(col.cols)
        col_path = mkpath(joinpath(path, string(col_name)))
        save_data(col, col_path)
    end
end

# =
# Delayed Collector
# = 

@kwdef struct DelayedCfg <: ColCfg 
    delay_time::Float64
end

mutable struct DelayedState{S} <: ColState 
    last_save_time::Float64
    sys_states::CircularBuffer{S}
end

struct DelayedCol{S} <: Collector 
    cfg::DelayedCfg
    state::DelayedState{S}
end
function get_collector(cfg::DelayedCfg, exp_cfg::ExperimentCfg, system::System, state=nothing)
    if isnothing(state)
        item = (deepcopy(system.time_info), deepcopy(system.state))
        past_states = CircularBuffer{typeof(item)}(2)
        push!(past_states, item)
        state = DelayedState(system.time_info.time, past_states)
    end

    DelayedCol(cfg, state)
end

function collect(col::DelayedCol, system::System)
    t = system.time_info.time
    col_state = col.state
    time_to_last_save = t - col_state.last_save_time
    if time_to_last_save > col.cfg.delay_time
        col_state.last_save_time = t
        push!(col_state.sys_states,
            (deepcopy(system.time_info), deepcopy(system.state))
        ) 
    end
end

function save_data(col::DelayedCol, path)
    times = [time_info.time for (time_info, _) in col.state.sys_states]
    idx = argmin(times)
    save_component_json(col.state.sys_states[idx][1], path, "time_info")
    save_component_serial(col.state.sys_states[idx][2], path, "state")
end

# =
# Checkpoint
# =

function check_checkpoint(cp::Nothing, experiment) end

function check_checkpoint(cfg::CheckpointCfg, experiment) 
    cp = experiment.checkpoint
    t = experiment.system.time_info.time
    if t - cp.last_time < cfg.delta_time
        return
    end
    cp.last_time = t

    cp_name = cp.last_used == "1" ? "2" : "1"
    cp_dir = joinpath(experiment.cfg.root, CHECKPOINT_DIRNAME) 
    path = joinpath(cp_dir, cp_name)

    cp.is_valid[cp_name] = false
    save_component_serial(experiment.system.state, path, "sys_state")
    serialize(joinpath(path, "sys_state/rng.bin"), experiment.system.rng)
    save_component_serial(experiment.system.time_info, path, "time_info")
    save_component_serial(experiment.col.state, path, "col_state")
    
    if cfg.save_data
        exp_cfg = experiment.cfg
        save_data(experiment.col, joinpath(exp_cfg.root, COL_DIRNAME))
    end

    cp.is_valid[cp_name] = true
    cp.last_used = cp_name

    serialize(joinpath(cp_dir, CHECKPOINT_INFO_NAME), cp)
end

function load_experiment(root, custom_step=nothing)
    cp_path = joinpath(root, CHECKPOINT_DIRNAME, CHECKPOINT_INFO_NAME)
    if !isfile(cp_path)
        error("Checkpoint info file does not exist at $cp_path")
    end
    cp_info::Checkpoint = deserialize(cp_path)

    cp_name = cp_info.last_used
    if !cp_info.is_valid[cp_name]
        cp_name = cp_name == "1" ? "2" : "1"
        if !cp_info.is_valid[cp_name]
            error("Every checkpoint is corrupted!")
        end
    end
        
    path = joinpath(root, CHECKPOINT_DIRNAME, cp_name)
    
    configs = convert(Dict{Symbol, Any}, JSON3.read(joinpath(root, EXP_COL_CONFIGS_NAME)))
    configs = load_dic_configs(configs)
    exp_cfg = configs[:experiment] 
    col_cfg = configs[:collector] 
    
    system = load_system(joinpath(root, "configs.json"), joinpath(path, "sys_state"), joinpath(path, "time_info"))
    col_state = load_component_serial(joinpath(path, "col_state"))

    col = get_collector(col_cfg, exp_cfg, system, col_state)

    return Experiment(exp_cfg, col, system, custom_step, cp_info)
end

include("rings/collectors.jl")

end # Collectors