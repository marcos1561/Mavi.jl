module Configs

export RingsCfg, has_types_cfg, get_spring_pars
export get_num_particles, get_area0, get_equilibrium_p0, get_particle_radius
export get_equilibrium_area, get_particles_area_contribution, get_ring_radius
export HarmTruncCfg
export InteractionMatrix, get_interaction_cfg, list_interactions, list_self_interactions
export IntCfg

using NLsolve

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

function MaviCfg.particle_radius(interaction_cfg::HarmTruncCfg) 
    return interaction_cfg.dist_eq / 2.0
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

@inline get_interaction_cfg(ring_id1, ring_id2, state, interaction::InteractionCfg) = interaction

@inline get_interaction_cfg(ring_type_1, ring_type_2, interaction::InteractionCfg) = interaction

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

@inline function get_ring_prop_by_name(dynamic_cfg::RingsCfg{U, T, I}, name, type=nothing) where {U<:Number, T, I}
    return getfield(dynamic_cfg, name)
end
@inline function get_ring_prop_by_name(dynamic_cfg::RingsCfg{U, T, I}, name, type) where {U<:AbstractVector, T, I}
    return getfield(dynamic_cfg, name)[type]
end

function get_rings_property(dynamic_cfg::RingsCfg, prop_func)
    prop_vec = Vector{Float64}(undef, dynamic_cfg.num_types)
    for i in 1:dynamic_cfg.num_types
        prop_vec[i] = prop_func(dynamic_cfg, i)
    end
    
    if length(prop_vec) == 1
        prop_vec = prop_vec[1]
    end

    return prop_vec
end

function get_spring_pars(dynamic_cfg::RingsCfg{U, T, F}, state=nothing, ring_id=nothing) where {U<:Number, T, F}
    return dynamic_cfg.k_spring, dynamic_cfg.l_spring
end

function get_spring_pars(dynamic_cfg::RingsCfg{U, T, F}, state, ring_id) where {U<:AbstractVector, T, F}
    t = state.types[ring_id]
    return dynamic_cfg.k_spring[t], dynamic_cfg.l_spring[t]
end


@inline function get_num_particles(dynamic_cfg::RingsCfg, type=nothing)
    return get_ring_prop_by_name(dynamic_cfg, :num_particles, type)
end

@inline function get_num_particles(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I}
    return get_rings_property(dynamic_cfg, get_num_particles)
end

@inline function get_particle_radius(dynamic_cfg::RingsCfg, type=nothing)
    inter = get_interaction_cfg(type, type, dynamic_cfg.interaction_finder)
    return particle_radius(inter)
end

function MaviCfg.particle_radius(dynamic_cfg::RingsCfg)
    return get_rings_property(dynamic_cfg, get_particle_radius)
    # for i in 1:dynamic_cfg.num_types
    #     inter = get_interaction_cfg(i, i, dynamic_cfg.interaction_finder)
    #     p_radius[i] = particle_radius(inter) 
    # end
    
    # if length(p_radius) == 1
    #     p_radius = p_radius[1]
    # end

    # return p_radius
end

"""
Returns the equilibrium area of the area force for
the given number of particles, taking into account `p0`.
"""
function get_area0(dynamic_cfg, type=nothing)
    num_particles = get_num_particles(dynamic_cfg, type)
    l_spring = get_ring_prop_by_name(dynamic_cfg, :l_spring, type)
    p0 = get_ring_prop_by_name(dynamic_cfg, :p0, type)
    return (num_particles * l_spring / p0)^2
end
get_area0(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_area0)

"""
Returns the p0 at which `area0` is equal to the equilibrium area
considering only the springs.
"""
function get_equilibrium_p0(dynamic_cfg::RingsCfg, type=nothing)
    num_particles = get_num_particles(dynamic_cfg, type)
    theta = 2 * π  / num_particles
    return 2 * (num_particles * (1 - cos(theta))/sin(theta))^.5
end
get_equilibrium_p0(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_equilibrium_p0)


"""
Equilibrium area of the polygon formed by the centers of the ring particles.

NOTE: The area contribution from the particles is not considered here.
"""
function get_equilibrium_area(dynamic_cfg::RingsCfg, type=nothing)
    l_spring = get_ring_prop_by_name(dynamic_cfg, :l_spring, type)

    num_particles = get_num_particles(dynamic_cfg, type) 
    a0 = get_area0(dynamic_cfg, type)

    p0_lim = get_equilibrium_p0(dynamic_cfg, type)
    a0_lim = (num_particles * l_spring / p0_lim)^2 

    if a0 < a0_lim
        return a0
    end
    theta = 2 * π / num_particles

    k_a = get_ring_prop_by_name(dynamic_cfg, :k_area, type)
    k_m = get_ring_prop_by_name(dynamic_cfg, :k_spring, type)

    function get_r(f)
        sqrt(f * 2 * a0 / (num_particles * sin(theta)))
    end

    function get_fm(f)
        k_m * (get_r(f) * sqrt(2 * (1 - cos(theta))) - l_spring)
    end

    function get_fm_total(f)
        2 * get_fm(f) * sin(theta / 2)
    end

    function get_fa(f)
        r = get_r(f)
        k_a * (a0 - 10/2 * r^2 * sin(theta)) * r * sin(theta)
    end

    function func!(F, f)
        F[1] = get_fa(f[1]) - get_fm_total(f[1])
    end

    sol = nlsolve(func!, [0.5])
    f_sol = sol.zero[1]

    return f_sol * a0
end
get_equilibrium_area(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_equilibrium_area)

"""
Contribution of the particles to the area of the rings.
The area of a ring is the area of the polygon formed by the centers of its particles (A_p),
plus the area of the particles that lies outside this polygon (A_c); this
function returns A_c.
"""
function get_particles_area_contribution(dynamic_cfg::RingsCfg, type=nothing)
    n = get_num_particles(dynamic_cfg, type)
    l_spring = get_ring_prop_by_name(dynamic_cfg, :l_spring, type)
    diameter = get_particle_radius(dynamic_cfg, type) * 2

    root_term = (diameter^2 - l_spring^2)^.5
    t = π/2 - atan(l_spring/root_term)
    area_intersect = 1/4 * (diameter^2 * t - l_spring * root_term)

    return n * π / 4 * diameter^2 * (1 - (n-2)/(2*n)) - n * area_intersect
end
get_particles_area_contribution(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_particles_area_contribution)

"""
Equilibrium radius of the ring, from the center to the outer edge of the particles.
It is calculated assuming its equilibrium area.
"""
function get_ring_radius(dynamic_cfg::RingsCfg, type=nothing)
    area_eq = get_equilibrium_area(dynamic_cfg, type)
    p_radius = get_particle_radius(dynamic_cfg, type)
    num_particles = get_num_particles(dynamic_cfg, type)
    radius_to_particle = (2 * area_eq / (num_particles * sin(2 * pi / num_particles)))^.5
    return radius_to_particle + p_radius
end
get_ring_radius(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_ring_radius)

end