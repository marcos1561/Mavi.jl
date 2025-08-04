module Rings

export RingsSystem, RingsState, Configs, InitStates, Integration
export get_num_particles, get_num_total_particles

using StaticArrays

import Mavi
import Mavi.Systems: get_num_total_particles

include("states.jl")
include("configs.jl")
include("neighbors.jl")
using .States
using .Configs
using .NeighborsMod

@kwdef struct RingsInfo{T, PN<:Union{ParticleNeighbors, Nothing}, RN<:Union{Neighbors, Nothing}, U}
    continuos_pos::Matrix{SVector{2, T}}
    areas::Vector{T}
    p_neigh::PN
    r_neigh::RN
    user_data::U
end

function RingsSystem(;state, space_cfg, dynamic_cfg, int_cfg, p_neighbors_cfg=nothing,
    r_neighbors=nothing, user_data=nothing)
    if has_types_cfg(dynamic_cfg) != has_types_func(state)
        if has_types_cfg(dynamic_cfg)
            error("DynamicCfg has multiple types, but state.types is nothing!")
        else
            error("DynamicCfg has only one type, but state.types is not nothing!")
        end
    end

    if !isnothing(r_neighbors) 
        num_max_neighbors = r_neighbors.only_count == false ? 15 : nothing 

        r_neighbors = Neighbors(
            num_entities=size(state.ring_pos, 3),
            num_max_neighbors=num_max_neighbors,
            device=int_cfg.device,
        )   
    end

    p_neighbors = nothing
    if !isnothing(p_neighbors_cfg) 
        num_max_neighbors = p_neighbors_cfg.only_count == false ? 15 : nothing 

        neigh = Neighbors(
            num_entities=length(state.pos),
            num_max_neighbors=num_max_neighbors,
            device=int_cfg.device,
            cfg=p_neighbors_cfg,
        )   
        p_neighbors = ParticleNeighbors(
            neighbors=neigh,
            type=p_neighbors_cfg.type,
            num_max_particles=num_max_particles(state),
        )
    end

    info = RingsInfo(
        continuos_pos=similar(state.rings_pos),
        areas=Vector{Float64}(undef, size(state.rings_pos, 2)),
        p_neigh=p_neighbors, r_neigh=r_neighbors,
        user_data=user_data,
    )

    Mavi.System(
        state=state, 
        space_cfg=space_cfg, 
        dynamic_cfg=dynamic_cfg, 
        int_cfg=int_cfg, 
        info=info, 
    )
end

@inline function get_num_particles(dynamic_cfg::RingsCfg{U, T, I}) where {U<:Number, T, I}
    return dynamic_cfg.num_particles
end

@inline function get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:AbstractVector, T, I}
    return dynamic_cfg.num_particles[state.types[ring_id]]
end

@inline function get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:Number, T, I}
    return get_num_particles(dynamic_cfg)
end

function Mavi.Systems.get_num_total_particles(system, state::RingsState)
    num_p = 0
    if state.types === nothing
        num_p = length(state.pos)
    else
        for t in state.types
            num_p += system.dynamic_cfg.num_particles[t]
        end
    end

    return num_p
end

function Mavi.Systems.particles_radius(system, dynamic_cfg::RingsCfg)
    state = system.state
    p_radius::Vector{Float64} = []
    for ring_id in get_active_ids(state)
        inter = get_interaction_cfg(ring_id, ring_id, state, dynamic_cfg.interaction_finder)
        num_p = get_num_particles(dynamic_cfg, state, ring_id)
        radius_i = Mavi.Configs.particle_radius(inter)
        for _ in 1:num_p
            push!(p_radius, radius_i)
        end
    end

    return p_radius
end

function Mavi.Systems.is_valid_pair(state::RingsState, dynamic_cfg::RingsCfg, i, j)
    num_max = num_max_particles(state)
    for idx in (i, j)
        p_id, ring_id = get_particle_ring_id(idx, num_max)
        num_ring_particles = get_num_particles(dynamic_cfg, state, ring_id)
        if p_id > num_ring_particles
            return false
        end
    end
    return true
end

include("integration.jl")
include("utils.jl")
include("init_states.jl")
# include("view.jl")
# include("space_checks.jl")

end