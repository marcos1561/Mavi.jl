module Mavi

export State, SecondLawState, SelfPropelledState, System

include("configs.jl")
using .Configs

abstract type State{T} end

"Particles state for Newton's second law (positions and velocities)"
@kwdef struct SecondLawState{T} <: State{T}
    pos::Matrix{T}
    vel::Matrix{T}
end

"Overdamped and self propelled state (positions and polarizations)"
@kwdef struct SelfPropelledState{T} <: State{T}
    pos::Matrix{T}
    pol_angle::Vector{T}
end

include("init_states.jl")
include("space_checks.jl")
include("chunks.jl")
using .SpaceChecks
using .ChunksMod

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
    Total force of all particles

    forces[d, i] = Total force on particle i in dimension d.
    """
    forces::Array{T, 2}
    
    """
    Distance between all particles

    dists[i, j] = Distance between particle i and j.
    """
    dists::Array{T, 2}

    "Number of particles"
    num_p::Int

    info::InfoT
    debug_info::DebugT
end
function System(;state::State{T}, space_cfg, dynamic_cfg, int_cfg, info=nothing, debug_info=nothing) where {T}
    all_inside, out_ids = check_inside(state, space_cfg.geometry_cfg)
    if all_inside == false
        throw("Particles with ids=$(out_ids) outside space.")
    end
    num_p = size(state.pos)[2]
    diffs = Array{T, 3}(undef, 2, num_p, num_p)
    forces = Array{T, 2}(undef, 2, num_p)
    dists = zeros(T, num_p, num_p)

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
            chunks_space_cfg, state, particle_radius(dynamic_cfg), extra_info=dynamic_cfg)
    end

    System(state, space_cfg, dynamic_cfg, int_cfg, chunks, diffs, forces, dists, num_p, info, debug_info)
end

include("integration.jl")
include("visualization.jl")
include("quantities.jl")

include("rings/rings.jl")

end
