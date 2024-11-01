module InitStates

using Mavi.Configs
using Mavi.Rings.Utils

function create_circle(center, radius, num_p)
    radius = radius
    pos = Matrix(undef, 2, num_p)

    theta = 2 * π / num_p
    for i in 1:num_p
        pos[1, i] = radius * cos((i - 1)* theta) + center[1]
        pos[2, i] = radius * sin((i - 1)* theta) + center[2]
    end

    return pos
end

function rectangular_grid(;
    num_cols, num_rows, num_particles, p_radius, 
    pad_x=0, pad_y=0, radius_k=1)
    ring_r = get_ring_radius(p_radius, num_particles)
    ring_length = 2 * (ring_r + p_radius)

    pad_x = pad_x * ring_r
    pad_y = pad_y * ring_r

    ring_pos = Array{Float64}(undef, 2, num_particles, num_cols * num_rows)
    idx = 1
    for col_id in 1:num_cols
        for row_id in 1:num_rows
            center_x = pad_x/2 + (col_id-1) * (ring_length + pad_x) + ring_length/2
            center_y = pad_y/2 + (row_id-1) * (ring_length + pad_x) + ring_length/2
            
            ring_pos[:, :, idx] = create_circle([center_x, center_y], ring_r * radius_k, num_particles)
            idx += 1
        end
    end

    space_l = num_cols * (pad_x + ring_length)
    space_h = num_rows * (pad_y + ring_length)
    geometry_cfg = RectangleCfg(space_l, space_h, (0, 0))

    return ring_pos, geometry_cfg
end

function random_pol(num_rings)
    return rand(Float64, num_rings) .* 2 * π
end

end