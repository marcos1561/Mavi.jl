module States

export RingsState
export get_active_ids

import Mavi

struct RingsState{T} <: Mavi.States.State{T}
    rings_pos::Array{T, 3} # size: (Nº of dimension, Nº of particles in the ring, Nº of rings)
    pos::Matrix{T} # size: (Nº of dimension , Total Nº of particles)
    pol::Vector{T}
    type::Vector{Int}
    active_ids::Vector{Int}
    uids::Vector{Int}
    num_active::Int
end
function RingsState(;rings_pos, pol, types)
    num_dims, num_max_particles, num_rings = size(rings_pos) 
    pos = reshape(rings_pos, num_dims, num_max_particles * num_rings)

    active_ids = Vector(1:num_rings)
    num_active = num_rings

    uids = Vector(1:num_rings)

    RingsState(rings_pos, pos, pol, types, active_ids, uids, num_active)
end

function get_active_ids(state::RingsState)
    return state.active_ids[1:state.num_active]
end

end