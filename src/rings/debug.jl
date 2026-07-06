module RingsDebug

export debug_pairwise_force!, debug_area_force!, debug_spring_force!
export DebugRingForces

using StaticArrays

using Mavi.Debug
using Mavi.Rings.States

abstract type RingsDebugInfo <: Debug.DebugInfo end

function debug_pairwise_force!(debug_info, f, ri, rj, i, j) end
function debug_area_force!(debug_info, f, ri, p_i) end
function debug_spring_force!(debug_info, f, ri, p_i1, p_i2) end

function debug_pairwise_force!(debug_info::ManyDebugs, args...)
    for d in debug_info.debugs
        debug_pairwise_force!(d, args...)
    end
end
function debug_area_force!(debug_info::ManyDebugs, args...) 
    for d in debug_info.debugs
        debug_area_force!(d, args...)
    end
end
function debug_spring_force!(debug_info::ManyDebugs, args...) 
    for d in debug_info.debugs
        debug_spring_force!(d, args...)
    end
end

# == 
# DebugRingForces
# ==

struct DebugRingForces{S, T} <: RingsDebugInfo
    ring_id::Int
    num_max_particles::Int
    area_forces::Vector{SVector{S, T}}
    spring_forces::Vector{SVector{S, T}}
    interaction_forces::Vector{Vector{SVector{S, T}}}
end
function DebugRingForces(state, ring_id)
    PosT = eltype(state.rings_pos) 
    np_max, nr = size(state.rings_pos)
    np = ring_num_particles(state, ring_id)

    interaction_forces = Vector{PosT}[]
    for _ in 1:nr
        push!(interaction_forces, Vector{PosT}(undef, np))
    end

    obj = DebugRingForces(ring_id, np_max, Vector{PosT}(undef, np), Vector{PosT}(undef, np), interaction_forces)
    clean_debug!(obj)
    return obj
end

function Debug.clean_debug!(debug_info::DebugRingForces)
    T = eltype(debug_info.area_forces)
    zero_el = zero(T) 
    debug_info.area_forces .= (zero_el,)
    debug_info.spring_forces .= (zero_el,)
    for forces_from_i in debug_info.interaction_forces
        forces_from_i .= (zero_el,)
    end
end

function debug_pairwise_force!(debug_info::DebugRingForces, f, ri, rj, i, j)
    if rj == debug_info.ring_id
        ri, rj = rj, ri
        i, j = j, i
        f = -f
    end

    if ri != debug_info.ring_id
        return
    end
    pid = get_particle_id(i, debug_info.num_max_particles, ri)
    debug_info.interaction_forces[rj][pid] += f
end

function debug_spring_force!(debug_info::DebugRingForces, f, ri, p_i1, p_i2)
    if ri != debug_info.ring_id
        return
    end
    debug_info.spring_forces[p_i1] += f
    debug_info.spring_forces[p_i2] -= f
end

function debug_area_force!(debug_info::DebugRingForces, f, ri, p_i)
    if ri != debug_info.ring_id
        return
    end
    debug_info.area_forces[p_i] += f
end

end # RingsDebug