module Systems

export System, particles_radius, get_forces, clean_forces!, get_num_total_particles, is_valid_pair, get_particle_radius
export get_particles_ids

using StaticArrays

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

"""
Base struct that represent a system of particles.

d=1: x axis
d=2: y axis
"""
struct System{T, ND, NT, StateT<:State{ND, T}, WallTypeT<:WallType, GeometryCfgT<:GeometryCfg, DynamicCfgT<:DynamicCfg, IntCfgT<:AbstractIntCfg,
    InfoT, DebugT}
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
    
    "Number of particles"
    num_p::Int

    time_info::TimeInfo

    info::InfoT
    debug_info::DebugT
end
function System(;state::State{ND, T}, space_cfg, dynamic_cfg, int_cfg, chunks=nothing, info=nothing, debug_info=nothing) where {ND, T}
    all_inside, out_ids = check_inside(state, space_cfg.geometry_cfg, get_particles_ids(state, info))
    if all_inside == false
        throw("Particles with ids=$(out_ids) outside space.")
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

    System(state, space_cfg, dynamic_cfg, int_cfg, chunks, forces, num_p, TimeInfo(0, 0.01), info, debug_info)
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

end