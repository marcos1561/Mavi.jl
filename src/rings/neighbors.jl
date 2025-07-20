module NeighborsMod

export NeighborsCfg, Neighbors, ParticleNeighbors 
export get_neigh, get_neigh_list, get_neigh_count, neigh_sum_buffers, neigh_clean!, neigh_update!

import Mavi

include("states.jl")
using .States

@kwdef struct NeighborsCfg 
    only_count::Bool = false
    type::Symbol = :all
end

abstract type AbstractNeighbors end

struct Neighbors{L<:Union{AbstractArray, Nothing}} <: AbstractNeighbors
    list::L
    count::Array{Int, 2}
    tol::Float64
end
function Neighbors(;num_entities, num_max_neighbors, device, tol=1.1)
    num_slices = 1
    if device isa Mavi.Configs.Threaded
        num_slices = Threads.nthreads()
    end

    neigh_count = fill(0, (num_entities, num_slices))
    neigh_array = nothing
    if !isnothing(num_max_neighbors)
        neigh_array = fill(-1, (num_max_neighbors, num_entities, num_slices))
    end

    Neighbors(neigh_array, neigh_count, tol)
end


"""
- type:  
    How neighbors are calculated.
    - `:rings` - Only particles from other rings are considered. Particles in 
    the same ring are not considered neighbors.
    - `:all` - Every particle is considered.
"""
@kwdef struct ParticleNeighbors{N<:Neighbors} <: AbstractNeighbors
    neighbors::N
    num_max_particles::Int
    type::Symbol = :rings
end

function get_neigh(neigh) end
get_neigh(neigh::Neighbors) = neigh
get_neigh(neigh::ParticleNeighbors) = neigh.neighbors

function get_neigh_list(neighbors, id) 
    neigh = get_neigh(neighbors)
    @view neigh.list[1:neigh.count[id, 1], id, 1]
end
get_neigh_count(neighbors) = @view get_neigh(neighbors).count[:, 1]

function neigh_clean!(neighbors) end

function neigh_clean!(neighbors::AbstractNeighbors)
   get_neigh(neighbors).count .= 0 
end

function neigh_sum_buffers(neighbors) end

function neigh_sum_buffers(neighbors::Neighbors{Nothing})
    neighbors.count[:, 1] .= sum(neighbors.count, dims=2)
end

function neigh_sum_buffers(neighbors::Neighbors{L}) where L<:AbstractArray
    neigh_main = @view neighbors.list[:, :, 1]
    count_main = @view neighbors.count[:, 1]

    for idx in 2:size(neighbors.list, 3)
        neigh = @view neighbors.list[:, :, idx]
        count = @view neighbors.count[:, idx]
        
        for ring_id in axes(neigh, 2)
            count_ring = count_main[ring_id]
            count_i = count[ring_id]
            for i in 1:count_i
                neigh_main[count_ring + i, ring_id] = neigh[i, ring_id] 
            end
            count_main[ring_id] += count_i
        end
    end
end

neigh_sum_buffers(neighbors::ParticleNeighbors) = neigh_sum_buffers(get_neigh(neighbors))

function neigh_update_data!(neighbors::Neighbors{Nothing}, i, j)
    count = neighbors.count
    tid = Threads.threadid()
    count[i, tid] += 1
    count[j, tid] += 1
end

function neigh_update_data!(neighbors::Neighbors, i, j)
    count = neighbors.count
    neigh_list = neighbors.list
    tid = Threads.threadid()

    count_i = count[i, tid] += 1
    count_j = count[j, tid] += 1
    
    neigh_list[count_i, i, tid] = j
    neigh_list[count_j, j, tid] = i
end

function neigh_update!(neighbors::Nothing, i, j, dist, max_dist) end

function neigh_update!(neighbors::Neighbors, i, j, dist, max_dist)
    if dist < max_dist * neighbors.tol
        neigh_update_data!(neighbors, i, j)    
    end
end

function neigh_update!(neighbors::Nothing, i, j, ri, rj, dist, max_dist) end

function neigh_update!(neighbors::ParticleNeighbors, i, j, ri, rj, dist, max_dist)
    inner_neigh = neighbors.neighbors
    t = neighbors.type
    # println("$(dist), $(max_dist), $(inner_neigh.tol)")

    if dist < max_dist * inner_neigh.tol && (t == :all || ri != rj)
        neigh_update_data!(inner_neigh, i, j)    
    end
end

end