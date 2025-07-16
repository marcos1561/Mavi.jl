module Rings

export RingsSystem, RingsState, Configs, InitStates, Integration
export get_num_particles

import Mavi

@kwdef struct RingsInfo{T}
    continuos_pos::Array{T, 3}
    areas::Vector{Float64}
end

@kwdef struct DebugInfo{T}
    graph_pos::Matrix{T}
    graph_radius::Vector{Float64}
    graph_color::Vector{Symbol}
    graph_type::Vector{Int}
end

include("states.jl")
include("configs.jl")
using .States
using .Configs

function RingsSystem(;state, space_cfg, dynamic_cfg, int_cfg)
    if has_types_cfg(dynamic_cfg) != has_types_func(state)
        if has_types_cfg(dynamic_cfg)
            error("DynamicCfg has multiple types, but state.types is nothing!")
        else
            error("DynamicCfg has only one type, but state.types is not nothing!")
        end
    end

    info = RingsInfo(
        continuos_pos = similar(state.rings_pos),
        areas = Vector{Float64}(undef, size(state.rings_pos)[3]),
    )

    num_total_p = length(state.pos)
    debug_info = DebugInfo(
        graph_pos = similar(state.pos),
        graph_radius = Vector{Float64}(undef, num_total_p),
        graph_color = Vector{Symbol}(undef, num_total_p),
        graph_type = Vector{Int}(undef, num_total_p),
    )

    Mavi.System(
        state=state, 
        space_cfg=space_cfg, 
        dynamic_cfg=dynamic_cfg, 
        int_cfg=int_cfg, 
        info=info, 
        debug_info=debug_info,
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