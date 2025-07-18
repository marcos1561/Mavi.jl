module Systems

export System, particles_radius, get_forces, clean_forces, get_num_total_particles

using Mavi.States: State
using Mavi.SpaceChecks
using Mavi.ChunksMod
using Mavi.Configs

mutable struct TimeInfo
    num_steps::Int
    time::Float64
end

"""
Base struct that represent a system of particles.

d=1: x axis
d=2: y axis
"""
struct System{T, StateT<:State{T}, WallTypeT<:WallType, GeometryCfgT<:GeometryCfg, DynamicCfgT<:DynamicCfg, IntCfgT<:AbstractIntCfg,
    InfoT, DebugT}
    state::StateT
    space_cfg::SpaceCfg{WallTypeT, GeometryCfgT}
    dynamic_cfg::DynamicCfgT
    int_cfg::IntCfgT
    chunks::Union{Chunks, Nothing}
    
    """
    Position difference between all particles
    
    diffs[d, i, j] = Position of particle i minus particle j in dimension d.

    d=1: x axis
    d=2: y axis
    """
    diffs::Array{T, 3}
    
    """
    Total force of all particles. One should use get_forces(system) to
    get a Num "Dimensions X Num Particles" matrix.

    forces[d, i, thread id] = Total force on particle i in dimension d.
    """
    forces::Array{T, 3}
    # forces_local::Union{Array{T, 3}, Nothing}
    
    """
    Distance between all particles

    dists[i, j] = Distance between particle i and j.
    """
    dists::Array{T, 2}

    "Number of particles"
    num_p::Int

    time_info::TimeInfo

    info::InfoT
    debug_info::DebugT
end
function System(;state::State{T}, space_cfg, dynamic_cfg, int_cfg, info=nothing, debug_info=nothing) where {T}
    all_inside, out_ids = check_inside(state, space_cfg.geometry_cfg)
    if all_inside == false
        throw("Particles with ids=$(out_ids) outside space.")
    end

    num_forces_slices = 1
    if int_cfg.device isa Threaded
        num_forces_slices = Threads.nthreads()
    end

    num_p = size(state.pos)[2]
    diffs = Array{T, 3}(undef, 2, num_p, num_p)
    forces = Array{T, 3}(undef, 2, num_p, num_forces_slices)
    dists = zeros(T, num_p, num_p)

    # forces_local = nothing
    # if int_cfg.device isa Threaded
    #     forces_local =  Array{T, 3}(undef, 2, num_p, num_forces_slices)
    # end

    chunks = nothing
    # if typeof(int_cfg) == ChunksIntCfg
    if has_chunks(int_cfg)
        chunks_cfg = int_cfg.chunks_cfg
        bounding_box = Configs.get_bounding_box(space_cfg.geometry_cfg) 
        chunks_space_cfg = SpaceCfg(
            wall_type=space_cfg.wall_type,
            geometry_cfg=bounding_box,
        )

        chunks = Chunks(chunks_cfg.num_cols, chunks_cfg.num_rows,
            chunks_space_cfg, state, minimum(particle_radius(dynamic_cfg)), extra_info=dynamic_cfg)
    end

    System(state, space_cfg, dynamic_cfg, int_cfg, chunks, diffs, forces, dists, num_p, TimeInfo(0, 0.01), info, debug_info)
end

@inline get_forces(system) = @view system.forces[:, :, 1]

@inline clean_forces(system) = system.forces .= 0

function particles_radius(system, dynamic_cfg)
    p_radius = particle_radius(dynamic_cfg)
    fill(p_radius, system.num_p)
end

@inline get_num_total_particles(system, state) = size(state.pos, 2)

@inline get_num_total_particles(system) = get_num_total_particles(system, system.state)
    
end