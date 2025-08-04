module States

export RingsState, get_active_ids, has_types_func
export num_max_particles, get_ring_id, get_particle_id, get_particle_ring_id, to_scalar_idx, get_inner_neigh_ids

using StaticArrays

import Mavi

mutable struct RingsState{T, VarT<:Union{Vector{Int}, Nothing}} <: Mavi.States.State{2, T}
    rings_pos::Matrix{SVector{2, T}} # size: (Nº of dimension, Nº of particles in the ring, Nº of rings)
    pos::Vector{SVector{2, T}} # size: (Nº of dimension , Total Nº of particles)
    pol::Vector{T}
    types::Union{Vector{Int}, Nothing}
    active_ids::VarT
    uids::VarT
    num_active::Int
end
function RingsState(;rings_pos, pol, types=nothing, is_variable=false)
    if eltype(rings_pos) <: Number
        num_dims, num_max_particles, num_rings = size(rings_pos) 
        rings_pos = copy(reinterpret(SVector{num_dims, eltype(rings_pos)}, vec(rings_pos)))
        rings_pos = reshape(rings_pos, (num_max_particles, num_rings))
    else
        num_max_particles, num_rings = size(rings_pos) 
    end

    pos = reshape(rings_pos, num_max_particles * num_rings)

    num_active = num_rings
    active_ids = nothing
    uids = nothing
    if is_variable
        active_ids = Vector(1:num_rings)
        uids = Vector(1:num_rings)
    end

    RingsState(rings_pos, pos, pol, types, active_ids, uids, num_active)
end

@inline get_active_ids(state::RingsState{U, Vector{Int}}) where U = state.active_ids[1:state.num_active]
@inline get_active_ids(state::RingsState{U, Nothing}) where U = axes(state.rings_pos, 2)

@inline has_types_func(state::RingsState) = !(state.types === nothing)

@inline is_variable_number(state::RingsState{U, Nothing}) where U = false
@inline is_variable_number(state::RingsState{U, Vector{Int}}) where U = true


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


end