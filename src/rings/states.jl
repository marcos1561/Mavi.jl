module States

export RingsState, ActiveState, ActiveRings, get_active_ids, has_types_func, add_ring
export num_max_particles, get_ring_id, get_particle_id, get_particle_ring_id, to_scalar_idx, get_inner_neigh_ids, calc_active_ids!
export get_num_active

using StaticArrays

import Mavi

struct ActiveState{M}
    mask::M
end
ActiveState() = ActiveState(true)

get_active_mask(active_state::ActiveState{Bool}, num_rings) = fill(active_state.mask, num_rings)
get_active_mask(active_state::ActiveState{T}, num_rings=nothing) where T <: AbstractVector = vec(active_state.mask)

@kwdef mutable struct ActiveRings
    mask::Vector{Bool}
    ids::Vector{Int}
    uids::Vector{Int}
    num_active::Int
    num_particles_active::Int
end

function calc_active_ids!(active) end
function calc_active_ids!(active::ActiveRings)
    pointer = 1
    ids = active.ids
    mask = active.mask
    for i in eachindex(active.mask)
        if mask[i]
            ids[pointer] = i
            pointer += 1 
        end
    end
    active.num_active = pointer - 1
end

struct RingsState{T, A<:Union{ActiveRings, Nothing}} <: Mavi.States.State{2, T}
    rings_pos::Matrix{SVector{2, T}} # size: (Nº of dimension, Nº of particles in the ring, Nº of rings)
    pos::Vector{SVector{2, T}} # size: (Nº of dimension , Total Nº of particles)
    pol::Vector{T}
    types::Union{Vector{Int}, Nothing}
    active::A
end
function RingsState(;rings_pos, pol, types=nothing, active_state=nothing)
    if eltype(rings_pos) <: Number
        num_dims, num_max_particles, num_rings = size(rings_pos) 
        rings_pos = copy(reinterpret(SVector{num_dims, eltype(rings_pos)}, vec(rings_pos)))
        rings_pos = reshape(rings_pos, (num_max_particles, num_rings))
    else
        num_max_particles, num_rings = size(rings_pos) 
    end

    pos = reshape(rings_pos, num_max_particles * num_rings)

    active_rings = nothing
    if !isnothing(active_state)
        active_rings = ActiveRings(
            mask=get_active_mask(active_state, num_rings), 
            ids=Vector(1:num_rings),
            uids=Vector(1:num_rings),
            num_active=0,
            num_particles_active=0,
        )
        calc_active_ids!(active_rings)
    end

    RingsState(rings_pos, pos, pol, types, active_rings)
end

@inline get_num_active(state::RingsState) = size(state.rings_pos, 2)
@inline get_num_active(state::RingsState{T, ActiveRings}) where T = state.active.num_active

@inline get_active_ids(state::RingsState) = axes(state.rings_pos, 2)
@inline get_active_ids(state::RingsState{T, ActiveRings}) where T = state.active.ids[1:state.active.num_active]

@inline has_types_func(state::RingsState) = !(state.types === nothing)

@inline is_variable_number(state::RingsState{U, Nothing}) where U = false
@inline is_variable_number(state::RingsState{U, ActiveRings}) where U = true


@inline num_max_particles(state::RingsState) = size(state.rings_pos, 1)

@inline function to_scalar_idx(state, ring_id, particle_id)
    (ring_id - 1) * num_max_particles(state) + particle_id
end

@inline get_ring_id(idx, num_max_particles) = div(idx-1, num_max_particles) + 1

@inline get_particle_id(idx, num_max_particles) = idx - (get_ring_id(idx, num_max_particles) - 1) * num_max_particles
@inline get_particle_id(idx, num_max_particles, ring_id) = idx - (ring_id - 1) * num_max_particles

function get_inner_neigh_ids(idx, num_max_particles)
    num_p = num_max_particles
    p_id = get_particle_id(idx, num_p)
    if p_id == 1 
        after_id = idx + 1
        before_id = idx + num_p - 1
    elseif p_id == num_p
        after_id = idx - num_p + 1
        before_id = idx - 1
    else
        after_id = idx + 1
        before_id = idx - 1
    end

    return before_id, after_id
end

@inline function get_particle_ring_id(idx, num_max_particles)
    ring_id = get_ring_id(idx, num_max_particles)
    particle_id = get_particle_id(idx, num_max_particles, ring_id)
    return particle_id, ring_id
end

function add_ring(state, pos, pol)
    active = state.active
    for i in eachindex(active.mask)
        if !active.mask[i]
            state.rings_pos[:, i] .= pos
            state.pol[i] = pol
            active.mask[i] = true
            active.uids[i] = maximum(active.uids) + 1
            active.num_active += 1
            return
        end
    end

    # @warn "Space not found to add a new Ring!"
end

end