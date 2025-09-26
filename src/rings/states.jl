module States

export RingsState, ActiveState
export FixRingsIds, VarRingsIds, get_rings_ids, get_num_active, get_particles_ids, get_num_total_particles, update_ids!
export has_types_func, add_ring!, remove_ring!
export num_max_particles, get_ring_id, get_particle_id, get_particle_ring_id, to_scalar_idx, get_inner_neigh_ids
export ring_num_particles, get_num_active

using StaticArrays

import Mavi.States as mv_states
using Mavi.States

# = 
# Rings Ids
# = 

abstract type RingsIds <: mv_states.AbstractParticleIds end

@kwdef mutable struct VarRingsIds <: RingsIds
    mask::Vector{Bool}
    ids::Vector{Int}
    uids::Vector{Int}
    p_ids::Vector{Int}
    num_active::Int
    num_p_active::Int
end

function mv_states.get_ids(rings_ids::VarRingsIds)
    ids = @view rings_ids.p_ids[1:rings_ids.num_p_active]
    return ids
end
mv_states.get_num(rings_ids::VarRingsIds) = rings_ids.num_p_active

function get_rings_ids(rings_ids::VarRingsIds)
    ids = @view rings_ids.ids[1:rings_ids.num_active]
    return ids
end
get_rings_num(rings_ids::VarRingsIds) = rings_ids.num_active 


@kwdef mutable struct FixRingsIds{RIDS, PIDS} <: RingsIds
    ids::RIDS
    p_ids::PIDS
end
function FixRingsIds(rings_pos, pos, types, num_particles)
    ids = axes(rings_pos, 2) 
    num_max_p = maximum(num_particles)
    if types isa AbstractVector
        p_ids = Int[]
        for ring_id in ids
            p_idx = to_scalar_idx(ring_id, 1, num_max_p) 
            for i in 1:ring_num_particles(num_particles, types, ring_id)
                push!(p_ids, p_idx + i - 1)
            end
        end
    else
        p_ids = eachindex(pos)
    end

    FixRingsIds(ids, p_ids)
end

mv_states.get_ids(rings_ids::FixRingsIds) = rings_ids.p_ids
mv_states.get_num(rings_ids::FixRingsIds) = length(rings_ids.p_ids)

get_rings_ids(rings_ids::FixRingsIds) = rings_ids.ids
get_rings_num(rings_ids::FixRingsIds) = length(rings_ids.ids)

# = 
# State
# = 

struct RingsState{T, RID<:Union{RingsIds, Nothing}, TP<:Union{Vector{Int}, Nothing}, NP<:Union{Vector{Int}, Int}} <: mv_states.State{2, T}
    rings_pos::Matrix{SVector{2, T}} # size: (Nº of dimension, Nº of particles in the ring, Nº of rings)
    pos::Vector{SVector{2, T}} # size: (Nº of dimension , Total Nº of particles)
    pol::Vector{T}
    num_particles::NP
    types::TP
    rings_ids::RID
end
function RingsState(;rings_pos, pol, num_particles=nothing, types=nothing, active_state=nothing)
    if num_particles isa AbstractVector && isnothing(types)
        if isnothing(types)
            error("argument 'types' is empty!")
        elseif length(num_particles) != length(types)
            ln = length(num_particles)
            lt = length(types)
            error("length(num_particles)=$ln must be equal to length(types)=$lt!")
        end
    end

    if eltype(rings_pos) <: Number
        num_dims, num_max_particles, num_rings = size(rings_pos) 
        rings_pos = copy(reinterpret(SVector{num_dims, eltype(rings_pos)}, vec(rings_pos)))
        rings_pos = reshape(rings_pos, (num_max_particles, num_rings))
    else
        num_max_particles, num_rings = size(rings_pos) 
    end

    pos = reshape(rings_pos, num_max_particles * num_rings)

    if isnothing(num_particles)
        num_particles = num_max_particles
    end

    rings_ids = nothing
    if !isnothing(active_state)
        rings_ids = VarRingsIds(
            mask=mv_states.get_active_mask(active_state, num_rings), 
            ids=Vector(1:num_rings),
            uids=Vector(1:num_rings),
            p_ids=Vector(1:num_rings*num_max_particles),
            num_active=0,
            num_p_active=0,
        )
    else
        rings_ids = FixRingsIds(rings_pos, pos, types, num_particles)
    end
    
    state = RingsState(rings_pos, pos, pol, num_particles, types, rings_ids)
    update_ids!(state)
    return state
end

@inline get_num_active(state::RingsState) = get_rings_num(state.rings_ids)
@inline get_rings_ids(state::RingsState) = get_rings_ids(state.rings_ids)
@inline get_rings_ids(system) = get_rings_ids(system.state.rings_ids)

@inline has_types_func(state::RingsState) = !(state.types === nothing)

@inline is_variable_number(state::RingsState{U, Nothing}) where U = false
@inline is_variable_number(state::RingsState{U, RingsIds}) where U = true

@inline num_max_particles(state::RingsState) = size(state.rings_pos, 1)

@inline function to_scalar_idx(ring_id, particle_id, num_particles)
    (ring_id - 1) * num_particles + particle_id
end
to_scalar_idx(state::RingsState, ring_id, particle_id) = to_scalar_idx(ring_id, particle_id, num_max_particles(state))

@inline get_ring_id(idx, num_max_particles) = div(idx-1, num_max_particles) + 1
@inline get_ring_id(state::RingsState, idx) = get_ring_id(idx, num_max_particles(state))

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

function add_ring!(state, pos, pol)
    rings_ids = state.rings_ids
    for i in eachindex(rings_ids.mask)
        if !rings_ids.mask[i]
            state.rings_pos[:, i] .= pos
            state.pol[i] = pol
            rings_ids.mask[i] = true
            rings_ids.uids[i] = maximum(rings_ids.uids) + 1
            rings_ids.num_active += 1
            return i
        end
    end
    return nothing
    # @warn "Space not found to add a new Ring!"
end

function remove_ring!(state, ring_id)
    rings_ids = state.rings_ids
    rings_ids.mask[ring_id] = false
    rings_ids.num_active -= 1
end

ring_num_particles(num_particles::Vector, ring_type) = num_particles[ring_type]
ring_num_particles(num_particles::Vector, types, rid) = num_particles[types[rid]]
ring_num_particles(num_particles::Int, types=nothing, rid=nothing) = num_particles
ring_num_particles(state::RingsState, ring_id) = ring_num_particles(state.num_particles, state.types, ring_id)

function calc_active_ids!(active, state) end 
function calc_active_ids!(active::VarRingsIds, state::RingsState)
    pointer = 1
    p_pointer = 1
    ids = active.ids
    p_ids = active.p_ids
    mask = active.mask
    for ring_id in eachindex(active.mask)
        if mask[ring_id]
            ids[pointer] = ring_id
            pointer += 1 
            
            ring_num_p = ring_num_particles(state, ring_id)
            scalar_id = to_scalar_idx(state, ring_id, 1)
            for pid in 1:ring_num_p
                p_ids[p_pointer] = scalar_id + pid - 1
                p_pointer += 1 
            end
        end
    end

    active.num_p_active = p_pointer - 1
    active.num_active = pointer - 1
end

function mv_states.update_ids!(state::RingsState)
    calc_active_ids!(state.rings_ids, state)
end

mv_states.get_particles_ids(state::RingsState) = mv_states.get_ids(state.rings_ids)
mv_states.get_num_total_particles(state::RingsState) = mv_states.get_num(state.rings_ids)

end