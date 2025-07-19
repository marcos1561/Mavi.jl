module Rings

export RingsSystem, RingsState, Configs, InitStates, Integration
export get_num_particles

import Mavi

include("states.jl")
include("configs.jl")
include("neighbors.jl")
using .States
using .Configs
using .NeighborsMod

@kwdef struct RingsInfo{T, PN<:Union{ParticleNeighbors, Nothing}, RN<:Union{Neighbors, Nothing}}
    continuos_pos::Array{T, 3}
    areas::Vector{T}
    p_neigh::PN
    r_neigh::RN
end

function RingsSystem(;state, space_cfg, dynamic_cfg, int_cfg, p_neighbors=nothing,
    r_neighbors=nothing,)
    if has_types_cfg(dynamic_cfg) != has_types_func(state)
        if has_types_cfg(dynamic_cfg)
            error("DynamicCfg has multiple types, but state.types is nothing!")
        else
            error("DynamicCfg has only one type, but state.types is not nothing!")
        end
    end

    if !isnothing(r_neighbors) 
        num_max_neighbors = r_neighbors.only_count == true ? 15 : nothing 

        r_neighbors = Neighbors(
            num_entities=size(state.ring_pos, 3),
            num_max_neighbors=num_max_neighbors,
            device=int_cfg.device,
        )   
    end

    if !isnothing(p_neighbors) 
        num_max_neighbors = p_neighbors.only_count == true ? 15 : nothing 

        neigh = Neighbors(
            num_entities=size(state.pos, 2),
            num_max_neighbors=num_max_neighbors,
            device=int_cfg.device,
        )   
        p_neighbors = ParticleNeighbors(
            neighbors=neigh,
            type=p_neighbors.type,
            num_max_particles=num_max_particles(state),
        )
    end

    info = RingsInfo(
        continuos_pos=similar(state.rings_pos),
        areas=Vector{Float64}(undef, size(state.rings_pos, 3)),
        p_neigh=p_neighbors, r_neigh=r_neighbors,
    )

    Mavi.System(
        state=state, 
        space_cfg=space_cfg, 
        dynamic_cfg=dynamic_cfg, 
        int_cfg=int_cfg, 
        info=info, 
    )
end

@inline function get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:AbstractVector, T, I}
    return dynamic_cfg.num_particles[state.types[ring_id]]
end

@inline function get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:Number, T, I}
    return dynamic_cfg.num_particles
end

@inline function get_num_particles(dynamic_cfg::RingsCfg{U, T, I}) where {U<:Number, T, I}
    return dynamic_cfg.num_particles
end

function Mavi.Systems.get_num_total_particles(system, state::RingsState)
    num_p = 0
    if state.types === nothing
        num_p = size(state.pos, 2)
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

include("integration.jl")
include("utils.jl")
include("init_states.jl")
include("view.jl")
# include("space_checks.jl")

end