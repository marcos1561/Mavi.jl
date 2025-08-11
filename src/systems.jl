module Systems

export System, particles_radius, get_forces, clean_forces!, get_num_total_particles, is_valid_pair, get_particle_radius
export get_particles_ids
export save, load_system

using StaticArrays, Serialization, JSON3, StructTypes

import Mavi
using Mavi.States
using Mavi.SpaceChecks
using Mavi.ChunksMod
using Mavi.Configs

function get_chunks(int_cfg::IntCfg, space_cfg::SpaceCfg, state, dynamic_cfg, extra_info=nothing)
    if !has_chunks(int_cfg)
        return nothing
    end

    chunks_cfg = int_cfg.chunks_cfg
    bounding_box = Configs.get_bounding_box(space_cfg.geometry_cfg) 
    chunks_space_cfg = SpaceCfg(
        wall_type=space_cfg.wall_type,
        geometry_cfg=bounding_box,
    )

    Chunks(chunks_cfg.num_cols, chunks_cfg.num_rows,
        chunks_space_cfg, state, minimum(Configs.particle_radius(dynamic_cfg)), extra_info=extra_info,
    )  
end

mutable struct TimeInfo
    num_steps::Int
    time::Float64
end

abstract type SystemType end
# StructTypes.StructType(::Type{S}) where S <: SystemType = StructTypes.Struct()
struct StandardSys <: SystemType end

"""
Base struct that represent a system of particles.

d=1: x axis
d=2: y axis
"""
struct System{T, ND, NT, StateT<:State{ND, T}, WallTypeT<:WallType, GeometryCfgT<:GeometryCfg, DynamicCfgT<:DynamicCfg, IntCfgT<:AbstractIntCfg,
    SysT<:SystemType, InfoT, DebugT}
    state::StateT
    space_cfg::SpaceCfg{WallTypeT, GeometryCfgT}
    dynamic_cfg::DynamicCfgT
    int_cfg::IntCfgT
    chunks::Union{Chunks, Nothing}
    
    """
    Total force of all particles. One should use get_forces(system) to
    get a Num "Dimensions X Num Particles" matrix.

    forces[d, i, thread id] = Total force on particle i in dimension d.
    """
    forces::SVector{NT, Vector{SVector{ND, T}}}
    # forces_local::Union{Array{T, 3}, Nothing}
    
    time_info::TimeInfo

    info::InfoT
    debug_info::DebugT
    type::SysT
end
function System(;state::State{ND, T}, space_cfg, dynamic_cfg, int_cfg, 
    chunks=nothing, info=nothing, debug_info=nothing, time_info=nothing, sys_type=:standard) where {ND, T}
    all_inside, out_ids = check_inside(state, space_cfg.geometry_cfg, get_particles_ids(state, info))
    if all_inside == false
        throw("Particles with ids=$(out_ids) outside space.")
    end

    if sys_type == :standard
        sys_type = StandardSys()
    end

    if isnothing(time_info)
        time_info = TimeInfo(0, 0.0)
    end

    num_forces_slices = 1
    if int_cfg.device isa Threaded
        num_forces_slices = Threads.nthreads()
    end

    # num_p = size(state.pos)[2]
    # forces = Array{T, 3}(undef, 2, num_p, num_forces_slices)
    
    num_p = length(state.pos)
    forces = SVector{num_forces_slices, Vector{SVector{ND, T}}}([
        Vector{SVector{ND, T}}(undef, num_p) for _ in 1:num_forces_slices
    ])

    # forces_local = nothing
    # if int_cfg.device isa Threaded
    #     forces_local =  Array{T, 3}(undef, 2, num_p, num_forces_slices)
    # end

    if isnothing(chunks)
        chunks = get_chunks(int_cfg, space_cfg, state, dynamic_cfg) 
    end

    if !isnothing(chunks)
        update_chunks!(chunks)
    end

    System(state, space_cfg, dynamic_cfg, int_cfg, chunks, forces, time_info, info, debug_info, sys_type)
end

# @inline get_forces(system) = @view system.forces[:, :, 1]
@inline get_forces(system) = system.forces[1]

@inline function clean_forces!(system)
    for f in system.forces
        f .= Scalar(zero(eltype(f)))
    end
end

function particles_radius(dynamic_cfg, state)
    p_radius = particle_radius(dynamic_cfg)
    fill(p_radius, get_num_total_particles(state))
end
particles_radius(system::System) = particles_radius(system.dynamic_cfg, system.state)

get_particle_radius(dynamic_cfg, state, idx) = particle_radius(dynamic_cfg)
get_particle_radius(system::System, idx) = get_particle_radius(system.dynamic_cfg, system.state, idx)

@inline get_num_total_particles(state::State, info=nothing) = length(state.pos)
@inline get_num_total_particles(system::System, state::State) = get_num_total_particles(state, system.info)
@inline get_num_total_particles(system::System) = get_num_total_particles(system.state, system.info)

@inline is_valid_pair(state::State, dynamic_cfg, i, j) = true
@inline is_valid_pair(system, i, j) = is_valid_pair(system.state, system.dynamic_cfg, i, j)

@inline get_particles_ids(state::State, dynamic_cfg=nothing) = eachindex(state.pos)
@inline get_particles_ids(system::System, state) = get_particles_ids(state, system.dynamic_cfg)
@inline get_particles_ids(system::System) = get_particles_ids(system.state, system.dynamic_cfg)



# ==
#  Save Load 
# ==

function save_component_serial(obj, root, name) 
    root_path = mkpath(joinpath(root, name))
    data = get_obj_save_data_json(obj)
    serialize(joinpath(root_path, "type.bin"), data[:type])
    serialize(joinpath(root_path, "data.bin"), data[:data])
end

# function save(obj, root, name)
#     serialize(joinpath(root, name * "_type.bin"), string(typeof(obj)))
#     save_component(obj, root, name)
# end

get_obj_save_data(obj) = obj
get_obj_save_data_json(obj) = (type=string(typeof(obj)), data=get_obj_save_data(obj))
# get_obj_type_and_data(obj::Tuple) = obj
get_obj_save_data_json(obj::Tuple) = [get_obj_save_data_json(x) for x in obj]

function save(system::System, root)
    mkpath(root)
    
    configs_data = Dict(
        "sys_type" => system.type,
        "space_cfg" => get_obj_save_data_json(system.space_cfg),
        "dynamic_cfg" => get_obj_save_data_json(system.dynamic_cfg),
        "int_cfg" => get_obj_save_data_json(system.int_cfg),
        "info" => get_obj_save_data_json(system.info),
        "debug_info" => get_obj_save_data_json(system.debug_info),
        "time_info" => get_obj_save_data_json(system.time_info),
    )

    open(joinpath(root, "configs.json"), "w") do io
        JSON3.pretty(io, configs_data)
    end
    save_component_serial(system.state, root, "state")
    # serialize(joinpath(root, "state.bin"), state_data)

    # serialize(joinpath(root, "sys_type.bin"), system.type)
    # save(system.state, root, "state")
    # save(system.space_cfg, root, "space_cfg")
    # save(system.dynamic_cfg, root, "dynamic_cfg")
    # save(system.int_cfg, root, "int_cfg")
    # save(system.info, root, "info")
    # save(system.debug_info, root, "debug_info")
    # save(system.time_info, root, "time_info")

    # serialize(joinpath(root, "pol.bin"), system.state.pol)
    # serialize(joinpath(root, "space_cfg.bin"), system.space_cfg)
    # serialize(joinpath(root, "dynamic_cfg.bin"), system.dynamic_cfg)
    # serialize(joinpath(root, "int_cfg.bin"), system.int_cfg)
end

get_serial_load_info(root) = (type=deserialize(joinpath(root, "type.bin")), data=root)

load_component(T::Type, data::JSON3.Object) = JSON3.read(JSON3.write(data), T)
load_component(T::Type, path) = deserialize(joinpath(path, "data.bin"))
function load_component(T::Type{Nothing}, data::JSON3.Array)
    data_loaded = []
    for value in data
        push!(data_loaded, load_component(get_saved_type(value), get_saved_data(value)))
    end
    return Tuple(data_loaded)
end

load_component(T::Type{Nothing}, _) = nothing

# function load_main_members(root)
#     return (
#         state = load_comp(root, "state"),
#         space_cfg = load_comp(root, "space_cfg"),
#         dynamic_cfg = load_comp(root, "dynamic_cfg"),
#         int_cfg = load_comp(root, "int_cfg"),
#         time_info = load_comp(root, "time_info"),
#     )
# end

get_saved_type(obj) = eval(Meta.parse(string(obj[:type])))
get_saved_type(obj::JSON3.Array) = Nothing

get_saved_data(obj) = obj[:data]
get_saved_data(obj::JSON3.Array) = obj

function load_dic_configs(configs)
    configs_loaded = Dict()
    for (name, value) in configs
        if name == :sys_type
            continue
        end
        println(name)
        println(typeof(value))
        configs_loaded[name] = load_component(get_saved_type(value), get_saved_data(value))
    end
    return configs_loaded
end

function load_system(configs, ::StandardSys)
    # main_members = load_main_members(root)
    # info = load_comp(root, "info")
    # debug_info = load_comp(root, "debug_info")

    configs_loaded = load_dic_configs(configs)
    # configs_loaded = Dict()
    # for (name, value) in configs
    #     if name == :sys_type
    #         continue
    #     end
        
    #     T = eval(Meta.parse(string(value[:type]))) 
    #     configs_loaded[name] = load_component(T, value[:data])
    # end

    # println(configs_loaded)

    System(
        state=configs_loaded[:state],
        space_cfg=configs_loaded[:space_cfg],
        dynamic_cfg=configs_loaded[:dynamic_cfg],
        int_cfg=configs_loaded[:int_cfg],
        info=configs_loaded[:info],
        time_info=configs_loaded[:time_info],
        debug_info=configs_loaded[:debug_info],
        # state=main_members.state,
        # space_cfg=main_members.space_cfg,
        # dynamic_cfg=main_members.dynamic_cfg,
        # int_cfg=main_members.int_cfg,
        # info=info,
        # time_info=main_members.time_info,
        # debug_info=debug_info,
    )
end

function load_system(root)
    configs = JSON3.read(joinpath(root, "configs.json"))
    configs = convert(Dict{Symbol, Any}, configs)
    configs[:state] = get_serial_load_info(joinpath(root, "state"))
    sys_type = eval(Meta.parse(configs[:sys_type]))
    # sys_type = deserialize(joinpath(root, "sys_type.bin"))
    
    load_system(configs, sys_type)
end

end