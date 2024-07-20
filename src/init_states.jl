module InitStates

using Mavi.Configs: RectangleCfg

"""
Return particles positions in a rectangle grid with spacing between
particles being `offset * radius`.

# Arguments
- num_p_x: Number of particle in x axis.
- num_p_y: Number of particle in y axis.
- offset: Space between particles, which is `offset * radius`.
- radius: Particle radius.

# Return
- pos::Matrix{Float64}
    Particle positions

- space_cfg::RectangleCfg
    SpaceCfg that contains all particles.
"""
function rectangular_grid(num_p_x, num_p_y, offset, radius)
    x = Vector{Float64}()
    y = Vector{Float64}()
    current_y = radius * (offset + 1)
    for i in 1:num_p_y
        current_x = -radius
        for j in 1:num_p_x
            push!(x, current_x + radius * (2 + offset))
            push!(y, current_y)
            current_x = x[end]
        end
        current_y += radius * (2 + offset)
    end
    
    space_cfg=RectangleCfg(
        length=num_p_x*2*radius + radius*offset*(num_p_x + 1),
        height=num_p_y*2*radius + radius*offset*(num_p_y + 1),
    )

    pos = Matrix{Float64}(undef, 2, num_p_x*num_p_y)
    pos[1, :] = x
    pos[2, :] = y
    return pos, space_cfg
end

"""
Create random velocitys for `num_p` particles with the maximum component 
value being `max_value`.
"""
function random_vel(num_p, max_value=1)
    return (rand(2, num_p).*2.0 .- 1.0) .* max_value
end

end