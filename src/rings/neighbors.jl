module NeighborsMod

export NeighborsCfg, Neighbors, ParticleNeighbors 
export get_neigh, get_neigh_list, get_neigh_count, neigh_sum_buffers, neigh_clean!, neigh_update!

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
    if device isa Mavi.ConfigsThreaded
        num_slices = Threads.nthreads()
    end

    neigh_count = fill(0, (num_entities, num_slices))
    neigh_array = nothing
    if !is_nothing(num_max_neighbors)
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

get_neigh_list(neighbors, id) = @view neighbors.neigh_array[1:neighbors.neigh_count[id, 1], id, 1]
get_neigh_count(neighbors, id) = neighbors.neigh_count[id, 1]

function neigh_clean!(neighbors) end

function neigh_clean!(neighbors::Neighbors)
   neighbors.count .= 0 
end

function neigh_sum_buffers(neighbors) end

function neigh_sum_buffers(neighbors::Neighbors{Nothing})
    neighbors.neigh_count[:, 1] .= sum(neighbors.neigh_array, dims=2)
end

function neigh_sum_buffers(neighbors::Neighbors{AbstractArray})
    neigh_main = @view neighbors.neigh_array[:, :, 1]
    count_main = @view neighbors.neigh_count[:, 1]

    for idx in 2:size(neighbors.neigh_array, 3)
        neigh = @view neighbors.neigh_array[:, :, idx]
        count = @view neighbors.neigh_count[:, idx]
        
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

neigh_adicional_check(neigh, i, j) = false

function neigh_adicional_check(neigh::ParticleNeighbors, i, j) 
    if neigh.type == :rings
        r1 = get_ring_id(i, neigh.num_max_particles)
        r2 = get_ring_id(j, neigh.num_max_particles)
        return r1 == r2
    else
        return false
    end
end


function neigh_update!(neighbors, inner_neigh, i, j, dist, max_dist) end

function neigh_update!(neighbors, inner_neigh::Neighbors{Nothing}, i, j, dist, max_dist)
    if dist < max_dist * inner_neigh.tol
        if neigh_adicional_check(neighbors, i, j)
            return
        end

        count = inner_neigh.count
        tid = Threads.threadid()
        count[i, tid] += 1
        count[j, tid] += 1
    end
end

function neigh_update!(neighbors, inner_neigh::Neighbors, i, j, dist, max_dist)
    if dist < max_dist * inner_neigh.tol
        if neigh_adicional_check(neighbors, i, j)
            return
        end

        count = inner_neigh.count
        neigh_list = inner_neigh.list
        tid = Threads.threadid()

        count_i = count[i, tid] += 1
        count_j = count[j, tid] += 1
        
        neigh_list[count_i, i, tid] = j
        neigh_list[count_j, j, tid] = i
    end
end

@inline function neigh_update!(neighbors, i, j, dist, max_dist)
    neigh_update!(neighbors, get_neigh(neighbors), i, j, dist, max_dist)    
end

end