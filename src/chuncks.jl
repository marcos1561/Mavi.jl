module ChuncksMod

export Chuncks, update_chuncks!

using Mavi: State
using Mavi.Configs: RectangleCfg

struct Chuncks{T}
    num_cols::Int
    num_rows::Int
    chunck_length::Float64
    chunck_height::Float64
    space_cfg::RectangleCfg
    steps_to_update::Int

    state::State{T}
    
    neighbors::Matrix{Vector{CartesianIndex{2}}}
    
    chunck_particles::Array{Int, 3}
    num_particles_in_chunck::Array{Int, 2}
end
function Chuncks(num_cols, num_rows, space_cfg::RectangleCfg, state, particle_r)
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
    
    chunck_length = space_cfg.length / num_cols
    chunck_height = space_cfg.height / num_rows

    nc = (ceil(0.5*chunck_length/particle_r) + 1) * ((ceil(0.5*chunck_height/particle_r) + 1))
    nc = trunc(Int, ceil(nc*1.1))
    chunck_particles = Array{Int}(undef, nc, num_rows, num_cols)
    num_particles_in_chunck = zeros(Int, num_rows, num_cols)

    Chuncks(num_cols, num_rows, chunck_length, chunck_height, space_cfg, 1, 
        state, neighbors, chunck_particles, num_particles_in_chunck)
end

"Update particles chunck positions."
function update_chuncks!(chuncks::Chuncks)
    state = chuncks.state

    space_h = chuncks.space_cfg.height
    chunck_l, chunck_h = chuncks.chunck_length, chuncks.chunck_height

    chuncks.num_particles_in_chunck .= 0

    num_p = size(state.pos)[2]
    for i in 1:num_p
        row_id = trunc(Int, div(-state.pos[2, i] + space_h, chunck_h)) + 1
        col_id = trunc(Int, div(state.pos[1, i], chunck_l)) + 1
        
        row_id -= row_id == (chuncks.num_rows + 1) ? 1 : 0
        col_id -= col_id == (chuncks.num_cols + 1) ? 1 : 0

        p_i = chuncks.num_particles_in_chunck[row_id, col_id] + 1
        chuncks.chunck_particles[p_i, row_id, col_id] = i 
        chuncks.num_particles_in_chunck[row_id, col_id] += 1
    end
end

end