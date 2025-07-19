module Configs

export RingsCfg, has_types_cfg
export HarmTruncCfg
export InteractionMatrix, get_interaction_cfg, list_interactions, list_self_interactions
export IntCfg

import Mavi.Configs as MaviCfg
# import Mavi.Configs: DynamicCfg, AbstractIntCfg, ChunksCfg, has_chunks, particle_radius

using Reexport
@reexport using Mavi.Configs


abstract type InteractionCfg end

@kwdef struct HarmTruncCfg <: InteractionCfg
    k_rep::Float64
    k_atr::Float64
    dist_eq::Float64
    dist_max::Float64
end

abstract type InteractionFinder{T} end

struct InteractionMatrix{T} <: InteractionFinder{T} 
    matrix::Matrix{T}
end

list_interactions(interactions::InteractionCfg) = [interactions]

function list_interactions(interactions::InteractionMatrix)
    mat = interactions.matrix
    return [mat[i, j] for i in 1:size(mat, 1) for j in i:size(mat, 2)]
end

list_self_interactions(interactions::InteractionCfg) = [interactions]

function list_self_interactions(interactions::InteractionMatrix)
    mat = interactions.matrix
    return [mat[i, i] for i in axes(mat, 1)]
end

@inline function get_interaction_cfg(ring_id1, ring_id2, state, interaction::InteractionMatrix)
    interaction.matrix[state.types[ring_id1], state.types[ring_id2]]
end

@inline function get_interaction_cfg(ring_type_1, ring_type_2, interaction::InteractionMatrix)
    interaction.matrix[ring_type_1, ring_type_2]
end

@inline function get_interaction_cfg(ring_id1, ring_id2, state, interaction::InteractionCfg) 
    interaction
end

@inline function get_interaction_cfg(ring_type_1, ring_type_2, interaction::InteractionCfg) 
    interaction
end

struct RingsCfg{U<:Union{AbstractVector, Number}, T<:InteractionCfg, InteracFinderT<:Union{InteractionFinder{T}, T}} <: DynamicCfg
    p0::U
    relax_time::U
    vo::U
    mobility::U
    rot_diff::U
    k_area::U
    k_spring::U
    l_spring::U
    num_particles::Union{Vector{Int}, Int}
    num_types::Int
    interaction_finder::InteracFinderT
end
function RingsCfg(;
    p0, relax_time, vo, mobility, rot_diff, k_area, k_spring, l_spring,
    num_particles, interaction_finder,
    )
    U = typeof(p0)
    if U <: Number && !(num_particles isa Int)
        throw(ArgumentError("If U is Number, num_particles must be Int"))
    elseif U <: AbstractVector && !(num_particles isa AbstractVector)
        throw(ArgumentError("If U is AbstractVector, num_particles must be Vector"))
    end

    if U <: Number
        num_types = 1
    else
        num_types = length(p0)
    end

    RingsCfg(
        p0, relax_time, vo, mobility, rot_diff, 
        k_area, k_spring, l_spring, num_particles, num_types, interaction_finder,
    )
end

@inline has_types_cfg(dynamic_cfg::RingsCfg{U, T, F}) where {U<:Number, T, F} = false
@inline has_types_cfg(dynamic_cfg::RingsCfg{U, T, F}) where {U<:AbstractVector, T, F} = true

function MaviCfg.particle_radius(dynamic_cfg::RingsCfg)
    p_radius = Vector{Float64}(undef, dynamic_cfg.num_types)
    for i in 1:dynamic_cfg.num_types
        inter = get_interaction_cfg(i, i, dynamic_cfg.interaction_finder)
        p_radius[i] = particle_radius(inter) 
    end
    
    if length(p_radius) == 1
        p_radius = p_radius[1]
    end

    return p_radius
end

function MaviCfg.particle_radius(interaction_cfg::HarmTruncCfg) 
    return interaction_cfg.dist_eq / 2.0
end

# @kwdef struct RingsIntCfg <: AbstractIntCfg
#     dt::Float64
#     chunks_cfg::Union{ChunksCfg, Nothing} = nothing
# end

# MaviCfg.has_chunks(int_cfg::RingsIntCfg) = !(int_cfg.chunks_cfg === nothing)

end