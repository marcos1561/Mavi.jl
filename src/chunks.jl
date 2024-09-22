module ChunksMod

export Chunks, update_chunks!

using Mavi: State
using Mavi.Configs: RectangleCfg

struct Chunks{T, StateT<:State{T}}
    num_cols::Int
    num_rows::Int
    chunk_length::Float64
    chunk_height::Float64
    space_cfg::RectangleCfg
    steps_to_update::Int

    state::StateT
    
    neighbors::Matrix{Vector{CartesianIndex{2}}}
    
    chunk_particles::Array{Int, 3}
    num_particles_in_chunk::Array{Int, 2}
end
function Chunks(num_cols, num_rows, space_cfg::RectangleCfg, state, particle_r)
    neighbors = Matrix{Vector{CartesianIndex{2}}}(undef, num_rows, num_cols)
    for i in 1:(num_rows-1)
        neighbors[i, 1] = [
            CartesianIndex(i+1, 1),
            CartesianIndex(i+1, 2),
            CartesianIndex(i, 2),
        ]
        
        for j in 2:(num_cols-1)
            neighbors[i, j] = [
                CartesianIndex(i+1, j-1),
                CartesianIndex(i+1, j),
                CartesianIndex(i+1, j+1),
                CartesianIndex(i, j+1),
            ]
        end

        neighbors[i, num_cols] = [
            CartesianIndex(i+1, num_cols-1),
            CartesianIndex(i+1, num_cols),
        ]
    end
    
    for j in 1:(num_cols-1)
        neighbors[num_rows, j] = [CartesianIndex(num_rows, j+1)]
    end
    neighbors[num_rows, num_cols] = []
    
    chunk_length = space_cfg.length / num_cols
    chunk_height = space_cfg.height / num_rows

    nc = (ceil(0.5*chunk_length/particle_r) + 1) * ((ceil(0.5*chunk_height/particle_r) + 1))
    nc = trunc(Int, ceil(nc*1.1))
    chunk_particles = Array{Int}(undef, nc, num_rows, num_cols)
    num_particles_in_chunk = zeros(Int, num_rows, num_cols)

    Chunks(num_cols, num_rows, chunk_length, chunk_height, space_cfg, 1, 
        state, neighbors, chunk_particles, num_particles_in_chunk)
end

"Update particles chunk positions."
function update_chunks!(chunks::Chunks)
    state = chunks.state

    space_h = chunks.space_cfg.height
    bottom_left = chunks.space_cfg.bottom_left
    chunk_l, chunk_h = chunks.chunk_length, chunks.chunk_height

    chunks.num_particles_in_chunk .= 0

    num_p = size(state.pos)[2]
    for i in 1:num_p
        row_id = trunc(Int, div(-state.pos[2, i] + bottom_left[2] + space_h, chunk_h)) + 1
        col_id = trunc(Int, div(state.pos[1, i] - bottom_left[1], chunk_l)) + 1
        
        row_id -= row_id == (chunks.num_rows + 1) ? 1 : 0
        col_id -= col_id == (chunks.num_cols + 1) ? 1 : 0

        p_i = chunks.num_particles_in_chunk[row_id, col_id] + 1
        chunks.chunk_particles[p_i, row_id, col_id] = i 
        chunks.num_particles_in_chunk[row_id, col_id] += 1
    end
end

end