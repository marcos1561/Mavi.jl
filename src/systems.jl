module Systems

export System, StandardSys
export particles_radius, get_forces, clean_forces!, get_num_total_particles, is_valid_pair, get_particle_radius, get_cm, reset_time

using StaticArrays, Serialization, JSON3, StructTypes, Random

import Mavi
using Mavi.States
using Mavi.SpaceChecks
using Mavi.ChunksMod
using Mavi.Configs
using Mavi.MovableObjects

function get_chunks(chunks_cfg::Union{ChunksCfg, Nothing}, space_cfg::SpaceCfg, pos, dynamic_cfg, extra_info=nothing)
    if isnothing(chunks_cfg)
        return nothing
    end

    bounding_box = Configs.get_bounding_box(space_cfg.geometry_cfg) 
    chunks_space_cfg = SpaceCfg(
        wall_type=space_cfg.wall_type,
        geometry_cfg=bounding_box,
    )

    Chunks(chunks_cfg.num_cols, chunks_cfg.num_rows,
        chunks_space_cfg, pos, minimum(Configs.particle_radius(dynamic_cfg)), extra_info=extra_info,
    )  
end

mutable struct TimeInfo
    num_steps::Int
    time::Float64
end
function reset(time_info::TimeInfo) 
    time_info.num_steps = 0
    time_info.time = 0
end

abstract type SystemType end
# StructTypes.StructType(::Type{S}) where S <: SystemType = StructTypes.Struct()
struct StandardSys <: SystemType end

"""
Base struct that represent a system of particles.

d=1: x axis
d=2: y axis
"""
struct System{T, ND, NT, 
    StateT<:State{ND, T}, WallTypeT<:WallType, GeometryCfgT<:GeometryCfg, DynamicCfgT<:DynamicCfg, IntCfgT<:AbstractIntCfg,
    ChunksT<:Union{Nothing, Chunks}, SysT<:SystemType, 
    SpaceDataT<:Union{Nothing, SpaceData}, InfoT, DebugT, RNGT<:AbstractRNG}
    state::StateT
    space_cfg::SpaceCfg{WallTypeT, GeometryCfgT}
    dynamic_cfg::DynamicCfgT
    int_cfg::IntCfgT
    chunks::ChunksT
    
    """
    Total force of all particles. One should use get_forces(system) to
    get the forces.

    forces[i, thread id] = Total force on particle i.
    """
    forces::SVector{NT, Vector{SVector{ND, T}}}
    # forces_local::Union{Array{T, 3}, Nothing}
    
    space_data::SpaceDataT

    time_info::TimeInfo

    info::InfoT
    debug_info::DebugT
    type::SysT
    rng::RNGT
end
function System(;state::State{ND, T}, space_cfg, dynamic_cfg, int_cfg, 
    chunks=nothing, info=nothing, debug_info=nothing, time_info=nothing, 
    sys_type=:standard, rng=nothing) where {ND, T}
    all_inside, out_ids = check_inside(state, space_cfg.geometry_cfg, get_particles_ids(state))
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

    num_p = length(state.pos)
    forces = SVector{num_forces_slices, Vector{SVector{ND, T}}}([
        Vector{SVector{ND, T}}(undef, num_p) for _ in 1:num_forces_slices
    ])

    if isnothing(chunks)
        chunks = get_chunks(int_cfg.chunks_cfg, space_cfg, state.pos, dynamic_cfg, state) 
    end

    if !isnothing(chunks)
        update_chunks!(chunks)
    end

    if isnothing(rng)
        rng=Random.GLOBAL_RNG
    end

    space_data = get_space_data(space_cfg)

    System(state, space_cfg, dynamic_cfg, int_cfg, chunks, forces, space_data, time_info, info, debug_info, sys_type, rng)
end

reset_time(system) = reset(system.time_info)

# @inline get_forces(system) = @view system.forces[:, :, 1]
@inline get_forces(system) = system.forces[1]

function clean_movable_objects_forces!(movable_objects) end
function clean_movable_objects_forces!(movable_objects::Vector{M}) where M <: MovableObject
    zero_force = zero(movable_objects[1].force)
    for object in movable_objects
        object.force = zero_force
    end
end
clean_movable_objects_forces!(system::System) = clean_movable_objects_forces!(get_movable_objects(system.state))

@inline function clean_forces!(system)
    for f in system.forces
        f .= Scalar(zero(eltype(f)))
    end
    clean_movable_objects_forces!(system)
end

function particles_radius(dynamic_cfg, state)
    p_radius = particle_radius(dynamic_cfg)
    if p_radius isa Number
        return fill(p_radius, get_num_total_particles(state))
    else
        return p_radius
    end
end
particles_radius(system::System) = particles_radius(system.dynamic_cfg, system.state)

get_particle_radius(dynamic_cfg, state, idx) = particle_radius(dynamic_cfg)
get_particle_radius(system::System, idx) = get_particle_radius(system.dynamic_cfg, system.state, idx)

@inline is_valid_pair(state::State, dynamic_cfg, i, j) = true
@inline is_valid_pair(system, i, j) = is_valid_pair(system.state, system.dynamic_cfg, i, j)

get_cm(system, state, entity_id) = state.pos[entity_id]

States.get_num_total_particles(system::System) = get_num_total_particles(system.state)

States.get_particles_ids(system::System) = get_particles_ids(system.state)
States.get_entities_ids(system::System) = get_entities_ids(system.state)

States.update_ids!(system::System) = update_ids!(system.state)

function system_deepcopy(system, sys_type)
    System(
        state=deepcopy(system.state),
        space_cfg=deepcopy(system.space_cfg),
        dynamic_cfg=deepcopy(system.dynamic_cfg),
        int_cfg=deepcopy(system.int_cfg),
        info=deepcopy(system.info),
        debug_info=deepcopy(system.debug_info),
        time_info=deepcopy(system.time_info),
        sys_type=sys_type,
        rng=deepcopy(system.rng),
    )
end

Base.deepcopy(system::System) = system_deepcopy(system, system.type)

end