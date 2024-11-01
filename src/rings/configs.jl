module Configs

export RingsCfg, particle_radius
export HarmTruncCfg
export get_interaction_cfg, InteractionMatrix
export RingsIntCfg, ChunksCfg

import Mavi.Configs: DynamicCfg, AbstractIntCfg, ChunksCfg, has_chunks, particle_radius

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

function get_interaction_cfg(t1, t2, interaction::InteractionMatrix)
    interaction.matrix[t1, t2]
end

struct RingsCfg{T, InteracFinderT<:InteractionFinder{T}} <: DynamicCfg
    p0::Vector{Float64}
    relax_time::Vector{Float64}
    vo::Vector{Float64}
    mobility::Vector{Float64}
    rot_diff::Vector{Float64}
    k_area::Vector{Float64}
    k_spring::Vector{Float64}
    l_spring::Vector{Float64}
    num_particles::Vector{Int}
    interaction_finder::InteracFinderT
    num_max_particles::Int
end
function RingsCfg(;
    p0,
    relax_time,
    vo,
    mobility,
    rot_diff,
    k_area,
    k_spring,
    l_spring,
    num_particles,
    interaction_finder,
    )
    num_max_particles = maximum(num_particles)
    RingsCfg(
        p0,
        relax_time,
        vo,
        mobility,
        rot_diff,
        k_area,
        k_spring,
        l_spring,
        num_particles,
        interaction_finder,
        num_max_particles,
    )
end

function particle_radius(dynamic_cfg::RingsCfg{HarmTruncCfg, I}) where I 
    return get_interaction_cfg(1, 1, dynamic_cfg.interaction_finder).dist_eq / 2.0
end

@kwdef struct RingsIntCfg <: AbstractIntCfg
    dt::Float64
    chunks_cfg::Union{ChunksCfg, Nothing} = nothing
end

has_chunks(int_cfg::RingsIntCfg) = !(int_cfg.chunks_cfg === nothing)

end