module States

export RingsState, get_active_ids, has_types_func
export num_max_particles, get_ring_id, get_particle_id, get_particle_ring_id, to_scalar_idx

import Mavi

mutable struct RingsState{T, VarT<:Union{Vector{Int}, Nothing}} <: Mavi.States.State{T}
    rings_pos::Array{T, 3} # size: (Nº of dimension, Nº of particles in the ring, Nº of rings)
    pos::Matrix{T} # size: (Nº of dimension , Total Nº of particles)
    pol::Vector{T}
    types::Union{Vector{Int}, Nothing}
    active_ids::VarT
    uids::VarT
    num_active::Int
end
function RingsState(;rings_pos, pol, types=nothing, is_variable=false)
    num_dims, num_max_particles, num_rings=size(rings_pos) 
    pos = reshape(rings_pos, num_dims, num_max_particles * num_rings)

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
@inline get_active_ids(state::RingsState{U, Nothing}) where U = axes(state.rings_pos, 3)

@inline has_types_func(state::RingsState) = !(state.types === nothing)

@inline is_variable_number(state::RingsState{U, Nothing}) where U = false
@inline is_variable_number(state::RingsState{U, Vector{Int}}) where U = true


@inline num_max_particles(state::RingsState) = size(state.rings_pos, 2)

@inline function to_scalar_idx(state, ring_id, particle_id)
    (ring_id - 1) * num_max_particles(state) + particle_id
end

@inline get_ring_id(idx, num_max_particles) = div(idx-1, num_max_particles) + 1

@inline get_particle_id(idx, num_max_particles) = idx - (get_ring_id(idx, num_max_particles) - 1) * num_max_particles
@inline get_particle_id(idx, num_max_particles, ring_id) = idx - (ring_id - 1) * num_max_particles

@inline function get_particle_ring_id(idx, num_particles)
    ring_id = get_ring_id(idx, num_particles)
    particle_id = get_particle_id(idx, num_particles, ring_id)
    return particle_id, ring_id
end


end