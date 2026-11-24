module NeighborsMod

export NeighborsCfg, HystericNeighborsCfg, Neighbors, RepeatableNeighbors, ParticleNeighbors
export get_neigh, get_neigh_list, get_neigh_count, neigh_sum_buffers, neigh_clean!, neigh_update!

import Mavi

include("states.jl")
using .States

@kwdef mutable struct NeighborsCfg 
    only_count::Bool = false
    type::Symbol = :all
    tol::Float64 = 1.1
end

abstract type AbstractNeighbors end

struct Neighbors{L<:Union{AbstractArray, Nothing}, D} <: AbstractNeighbors
    list::L
    count::Array{Int, 2}
    cfg::NeighborsCfg
    device::D
end
function Neighbors(;num_entities, num_max_neighbors, device, cfg)
    num_slices = 1
    if device isa Mavi.Configs.Threaded
        num_slices = Threads.nthreads()
    end

    neigh_count = fill(0, (num_entities, num_slices))
    neigh_array = nothing
    if !isnothing(num_max_neighbors)
        neigh_array = fill(-1, (num_max_neighbors, num_entities, num_slices))
    end

    Neighbors(neigh_array, neigh_count, cfg, device)
end

get_slide_idx(::Mavi.Configs.Sequencial) = 1
get_slide_idx(::Mavi.Configs.Threaded) = Mavi.Integration.get_thread_index()

"""
    RepeatableNeighbors(neighbors::Neighbors)
    RepeatableNeighbors(; kwargs...)

Accept repeated observations of undirected pairs and expose unique neighbors
after `neigh_sum_buffers`. Self-pairs are ignored. Each worker records pairs in
its own set; reduction unions those sets before building lists and counts.
Call `neigh_clean!` before recording a new step. Lists grow when necessary.
Keyword arguments are forwarded to `Neighbors`.
"""
mutable struct RepeatableNeighbors{N<:Neighbors} <: AbstractNeighbors
    neighbors::N
    pair_buffers::Vector{Set{Tuple{Int, Int}}}
    unique_pairs::Set{Tuple{Int, Int}}
end
function RepeatableNeighbors(neighbors::Neighbors)
    buffers = [Set{Tuple{Int, Int}}() for _ in axes(neighbors.count, 2)]
    RepeatableNeighbors(neighbors, buffers, Set{Tuple{Int, Int}}())
end
RepeatableNeighbors(; kwargs...) = RepeatableNeighbors(Neighbors(; kwargs...))

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

@kwdef struct HystericNeighborsCfg
    neigh_cfg::NeighborsCfg
    offset_k::Float64
end

struct HystericNeighbors{N} <: AbstractNeighbors
    cfg::HystericNeighborsCfg
    neighbors::N
    is_neigh_buffer::BitMatrix
end
HystericNeighbors(cfg, neighbors, num_entities) = HystericNeighbors(cfg, neighbors, falses(num_entities, num_entities))

function get_neigh(neigh) end
get_neigh(neigh::Neighbors) = neigh
get_neigh(neigh::Union{ParticleNeighbors, HystericNeighbors, RepeatableNeighbors}) = neigh.neighbors

function get_neigh_list(neighbors, id) 
    neigh = get_neigh(neighbors)
    @view neigh.list[1:neigh.count[id, 1], id, 1]
end
get_neigh_count(neighbors) = @view get_neigh(neighbors).count[:, 1]

function neigh_clean!(neighbors) end

function neigh_clean!(neighbors::AbstractNeighbors)
    get_neigh(neighbors).count .= 0 
end

function neigh_clean!(neighbors::RepeatableNeighbors)
    neighbors.neighbors.count .= 0
    for pairs in neighbors.pair_buffers
        empty!(pairs)
    end
    empty!(neighbors.unique_pairs)
    return nothing
end

function neigh_sum_buffers(neighbors) end

function neigh_sum_buffers(neighbors::Neighbors{Nothing, D}) where D
    neighbors.count[:, 1] .= sum(neighbors.count, dims=2)
end

function neigh_sum_buffers(neighbors::Neighbors{L, D}) where {L<:AbstractArray, D}
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

neigh_sum_buffers(neighbors::AbstractNeighbors) = neigh_sum_buffers(get_neigh(neighbors))

function neigh_sum_buffers(neighbors::RepeatableNeighbors)
    pairs = neighbors.unique_pairs
    empty!(pairs)
    for buffer in neighbors.pair_buffers
        union!(pairs, buffer)
    end

    inner = neighbors.neighbors
    inner.count .= 0
    for (i, j) in pairs
        inner.count[i, 1] += 1
        inner.count[j, 1] += 1
    end

    # Dynamically changes the size of the neighbors list buffer 
    # if !isnothing(inner.list)
    #     required = maximum(@view(inner.count[:, 1]); init=0)
    #     if required > size(inner.list, 1)
    #         capacity = max(required, 2 * size(inner.list, 1))
    #         list = fill(-1, capacity, size(inner.list, 2), size(inner.list, 3))
    #         inner = Neighbors(list, inner.count, inner.cfg, inner.device)
    #         neighbors.neighbors = inner
    #     end

    #     inner.count .= 0
    #     for (i, j) in pairs
    #         count_i = inner.count[i, 1] += 1
    #         count_j = inner.count[j, 1] += 1
    #         inner.list[count_i, i, 1] = j
    #         inner.list[count_j, j, 1] = i
    #     end
    #     for i in axes(inner.count, 1)
    #         sort!(@view inner.list[1:inner.count[i, 1], i, 1])
    #     end
    # end
    
    return nothing
end

function neigh_update_data!(neighbors::Neighbors{Nothing, D}, i, j) where D
    count = neighbors.count
    idx = get_slide_idx(neighbors.device)
    count[i, idx] += 1
    count[j, idx] += 1
end

function neigh_update_data!(neighbors::Neighbors, i, j)
    count = neighbors.count
    neigh_list = neighbors.list
    idx = get_slide_idx(neighbors.device)

    count_i = count[i, idx] += 1
    count_j = count[j, idx] += 1
    
    neigh_list[count_i, i, idx] = j
    neigh_list[count_j, j, idx] = i
end

function neigh_update!(neighbors::Nothing, i, j, dist, max_dist) end

function neigh_update!(neighbors::Neighbors, i, j, dist, max_dist)
    if dist < max_dist * neighbors.cfg.tol
        neigh_update_data!(neighbors, i, j)    
    end
end

function neigh_update!(neighbors::RepeatableNeighbors, i, j, dist, max_dist)
    inner = neighbors.neighbors
    if i != j && dist < max_dist * inner.cfg.tol
        idx = get_slide_idx(inner.device)
        push!(neighbors.pair_buffers[idx], minmax(i, j))
    end
    return nothing
end

function neigh_update!(neighbors::HystericNeighbors, i, j, dist, max_dist)
    inner_neigh = neighbors.neighbors
    is_neigh_buffer = neighbors.is_neigh_buffer
    offset_k = neighbors.cfg.offset_k
    max_dist = inner_neigh.cfg.tol * max_dist
    if is_neigh_buffer[i, j]
        if dist < max_dist * (1 + offset_k)
            neigh_update_data!(inner_neigh, i, j)
        else
            is_neigh_buffer[i, j] = false
            is_neigh_buffer[j, i] = false
        end
    else
        if dist < max_dist
            is_neigh_buffer[i, j] = true
            is_neigh_buffer[j, i] = true
            neigh_update_data!(inner_neigh, i, j)
        end
    end
end

function neigh_update!(neighbors::Nothing, i, j, ri, rj, dist, max_dist) end

function neigh_update!(neighbors::ParticleNeighbors, i, j, ri, rj, dist, max_dist)
    inner_neigh = neighbors.neighbors
    t = neighbors.type
    # println("$(dist), $(max_dist), $(inner_neigh.tol)")

    if dist < max_dist * inner_neigh.cfg.tol && (t == :all || ri != rj)
        neigh_update_data!(inner_neigh, i, j)    
    end
end

end
