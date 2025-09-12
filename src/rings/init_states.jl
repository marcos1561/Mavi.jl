module InitStates

export rectangular_grid, random_pol, create_circle

using StaticArrays

using Mavi.Configs
using Mavi.Rings.Utils

function create_circle(center, radius, num_p; dtype=Float64, to_svector=false)
    radius = radius
    pos = Matrix{dtype}(undef, 2, num_p)

    theta = 2 * π / num_p
    for i in 1:num_p
        pos[1, i] = radius * cos((i - 1)* theta) + center[1]
        pos[2, i] = radius * sin((i - 1)* theta) + center[2]
    end

    if to_svector
        pos = copy(reinterpret(SVector{2, eltype(pos)}, vec(pos)))
    end

    return pos
end

function rectangular_grid(;
    num_cols, num_rows, num_particles, p_radius, types=nothing, 
    pad_x=0, pad_y=0, radius_k=1)
    if types === nothing
        num_particles = [num_particles]
        p_radius = [p_radius]
        types = fill(1, num_cols * num_rows)
    end

    ring_r = get_ring_radius.(p_radius, num_particles)
    ring_length = 2 * (ring_r + p_radius)

    max_ring_r = maximum(ring_r)
    max_ring_length = maximum(ring_length)
    max_num_particles = maximum(num_particles)

    pad_x = pad_x * max_ring_r
    pad_y = pad_y * max_ring_r

    ring_pos = Array{Float64}(undef, 2, max_num_particles, num_cols * num_rows)
    ring_pos .= 0
    idx = 1
    for col_id in 1:num_cols
        for row_id in 1:num_rows
            center_x = pad_x/2 + (col_id-1) * (max_ring_length + pad_x) + max_ring_length/2
            center_y = pad_y/2 + (row_id-1) * (max_ring_length + pad_x) + max_ring_length/2
            
            ring_type = types[idx]
            num_p = num_particles[ring_type]
            ring_pos[:, 1:num_p, idx] = create_circle([center_x, center_y], ring_r[ring_type] * radius_k, num_p)
            idx += 1
        end
    end

    space_l = num_cols * (pad_x + max_ring_length)
    space_h = num_rows * (pad_y + max_ring_length)
    geometry_cfg = RectangleCfg(length=space_l, height=space_h)

    return ring_pos, geometry_cfg
end

function random_pol(num_rings)
    return rand(Float64, num_rings) .* 2 * π
end

end