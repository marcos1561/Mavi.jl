module Collectors

export ExperimentCfg, Experiment, run_experiment
export DelayedCfg

using Serialization

using Mavi.Systems
using Mavi.MaviSerder
using Mavi.Rings.States
using Mavi.Rings.Integration

using Mavi.Utils.Progress

# Module Configs
COL_DIRNAME = "data"
FINAL_STATE_NAME = "final_state"

abstract type ColCfg end
abstract type ColState end
abstract type Collector end

function get_collector(col_cfg::ColCfg, exp_cfg, system) end
function collect(col::Collector, system) end
function save_data(col::Collector, path) end

# = 
# Experiment
# = 

@kwdef struct ExperimentCfg 
    tf::Float64
    root::String
    save_final_state::Bool=true
end

struct Experiment{S<:System, F<:Union{Nothing, Function}} 
    cfg::ExperimentCfg
    col::Collector
    system::S
    custom_step::F
end
function Experiment(cfg::ExperimentCfg, col_cfg::ColCfg, system::System, custom_step=nothing)
    Experiment(cfg, get_collector(col_cfg, cfg, system), system, custom_step)
end

function run_experiment(experiment::Experiment)
    system = experiment.system
    col = experiment.col
    cfg = experiment.cfg

    save_system_configs(system, cfg.root)

    col_path = mkpath(joinpath(cfg.root, COL_DIRNAME))

    experiment_step! = experiment.custom_step
    if isnothing(experiment_step!)
        experiment_step! = step!
    end

    prog = ProgContinuos(init=system.time_info.time, final=cfg.tf)
    while system.time_info.time < cfg.tf
        experiment_step!(system)
        collect(col, system)
        show_progress(prog, system.time_info.time)
    end
    save_data(col, col_path)
    
    if cfg.save_final_state
        save_system(system, joinpath(cfg.root, FINAL_STATE_NAME))
        # save_component_serial(system.state, cfg.root, FINAL_STATE_NAME)
    end

    show_finish(prog)
end

# =
# Delayed Collector
# = 

@kwdef struct DelayedCfg <: ColCfg 
    delay_time::Float64
end

mutable struct DelayedState <: ColState 
    last_save_time::Float64
end

struct DelayedCol <: Collector 
    cfg::DelayedCfg
    state::DelayedState
    path::String
end
function get_collector(cfg::DelayedCfg, exp_cfg::ExperimentCfg, system::System)
    DelayedCol(cfg, DelayedState(system.time_info.time), joinpath(exp_cfg.root, COL_DIRNAME))
end

function collect(col::DelayedCol, system::System)
    t = system.time_info.time
    col_state = col.state
    time_to_last_save = t - col_state.last_save_time
    if time_to_last_save > col.cfg.delay_time
        col_state.last_save_time = t
        save_component_serial(system.state, col.path, "state")
    end
end

end # Collectors