module ChunksMod

export Chunks, update_chunks!, get_chunk_particles, get_chunk_rect

using StaticArrays

using Mavi.States: State
using Mavi.Configs

struct Chunks{N, T, P, InfoT}
    num_cols::Int
    num_rows::Int
    chunk_length::Float64
    chunk_height::Float64
    geometry_cfg::RectangleCfg{N, T}
    steps_to_update::Int

    pos::P
    extra_info::InfoT

    neighbors::Matrix{Vector{CartesianIndex{2}}}
    
    chunk_particles::Array{Int, 3}
    num_particles_in_chunk::Array{Int, 2}
end
function Chunks(num_cols, num_rows, space_cfg::SpaceCfg{W, RectangleCfg{N, T}}, pos, particle_r; extra_info=nothing) where {W, N, T}
    neighbors = get_neighbors(num_rows, num_cols, get_main_wall(space_cfg.wall_type))
    
    chunk_length = space_cfg.geometry_cfg.length / num_cols
    chunk_height = space_cfg.geometry_cfg.height / num_rows

    nc = (ceil(0.5*chunk_length/particle_r) + 1) * ((ceil(0.5*chunk_height/particle_r) + 1))
    # nc = trunc(Int, ceil(nc*1.1))
    nc = trunc(Int, ceil(nc*2))
    chunk_particles = Array{Int}(undef, nc, num_rows, num_cols)
    num_particles_in_chunk = zeros(Int, num_rows, num_cols)

    Chunks(num_cols, num_rows, chunk_length, chunk_height, space_cfg.geometry_cfg, 1, 
        pos, extra_info, neighbors, chunk_particles, num_particles_in_chunk)
end

get_chunk_particles(chunks::Chunks, id) = get_chunk_particles(chunks, Tuple(id)...)
function get_chunk_particles(chunks::Chunks, row_id, col_id)
    num_p = chunks.num_particles_in_chunk[row_id, col_id]
    pids = @view chunks.chunk_particles[1:num_p, row_id, col_id]
    return pids
end

get_chunk_rect(chunks::Chunks, id) = get_chunk_rect(chunks, Tuple(id)...)
function get_chunk_rect(chunks::Chunks, row_id, col_id)
    up, right = SVector(0, 1), SVector(1, 0) 
    tl = chunks.geometry_cfg.bottom_left + chunks.geometry_cfg.height * up
    l, h = chunks.chunk_length, chunks.chunk_height
    return RectangleCfg(
        length=l,
        height=h,
        bottom_left=tl - row_id * h * up + (col_id - 1) * l * right 
    )
end

function get_neighbors(num_rows, num_cols, wall_type::PeriodicWalls)
    neighbors = Matrix{Vector{CartesianIndex{2}}}(undef, num_rows, num_cols)

    function get_id(x, num_t)
        if x == 0
            return num_t
        end

        if x % (num_t + 1) == 0
            return 1
        end
        return x
    end

    for i in 1:num_rows
        for j in 1:num_cols
            neighbors[i, j] = [
                CartesianIndex(get_id(i+1, num_rows), get_id(j  , num_cols)),
                CartesianIndex(get_id(i+1, num_rows), get_id(j+1, num_cols)),
                CartesianIndex(get_id(i  , num_rows), get_id(j+1, num_cols)),
                CartesianIndex(get_id(i-1, num_rows), get_id(j+1, num_cols)),
            ]
        end
    end

    return neighbors
end

function get_neighbors(num_rows, num_cols, wall_type::Union{RigidWalls, SlipperyWalls})
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
    return neighbors
end

function update_particle_chunk!(chunks, i)
    pos = chunks.pos

    space_h = chunks.geometry_cfg.height
    bottom_left = chunks.geometry_cfg.bottom_left
    chunk_l, chunk_h = chunks.chunk_length, chunks.chunk_height

    
    pos_i = pos[i]
    row_id = trunc(Int, div(-pos_i[2] + bottom_left[2] + space_h, chunk_h)) + 1
    col_id = trunc(Int, div(pos_i[1] - bottom_left[1], chunk_l)) + 1
    # row_id = trunc(Int, div(-pos[2, i] + bottom_left[2] + space_h, chunk_h)) + 1
    # col_id = trunc(Int, div(pos[1, i] - bottom_left[1], chunk_l)) + 1
    
    # if occursin("RingsState", string(typeof(chunks.extra_info)))
    #     println(pos_i)
    #     println(row_id, ", ", col_id)
    #     println("====")
    # end

    row_id -= row_id == (chunks.num_rows + 1) ? 1 : 0
    col_id -= col_id == (chunks.num_cols + 1) ? 1 : 0

    p_i = chunks.num_particles_in_chunk[row_id, col_id] + 1
    chunks.chunk_particles[p_i, row_id, col_id] = i 
    chunks.num_particles_in_chunk[row_id, col_id] += 1
end

"Update particles chunk positions."
function update_chunks!(chunks::Chunks)
    # num_p = size(chunks.pos)[2]
    num_p = length(chunks.pos)
    chunks.num_particles_in_chunk .= 0
    for i in 1:num_p
        update_particle_chunk!(chunks, i)
    end
end

function update_chunks!(chunks::Nothing) end

# function update_chunks!(chunks::Chunks)
#     pos = chunks.pos

#     space_h = chunks.geometry_cfg.height
#     bottom_left = chunks.geometry_cfg.bottom_left
#     chunk_l, chunk_h = chunks.chunk_length, chunks.chunk_height

#     chunks.num_particles_in_chunk .= 0

#     num_p = size(pos)[2]
#     for i in 1:num_p
#         row_id = trunc(Int, div(-pos[2, i] + bottom_left[2] + space_h, chunk_h)) + 1
#         col_id = trunc(Int, div(pos[1, i] - bottom_left[1], chunk_l)) + 1
        
#         row_id -= row_id == (chunks.num_rows + 1) ? 1 : 0
#         col_id -= col_id == (chunks.num_cols + 1) ? 1 : 0

#         p_i = chunks.num_particles_in_chunk[row_id, col_id] + 1
#         chunks.chunk_particles[p_i, row_id, col_id] = i 
#         chunks.num_particles_in_chunk[row_id, col_id] += 1
#     end
# end

end