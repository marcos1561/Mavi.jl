module InitStates

export rectangular_grid, random_vel

using StaticArrays
using Random

using Mavi.Configs: RectangleCfg

"""
Return particles positions in a rectangle grid with spacing between
particles being `offset * radius`.

# Arguments
- num_p_x:  
    Number of particle in x axis.

- num_p_y:   
    Number of particle in y axis.

- offset:   
    Space between particles, which is `offset * radius`.

- radius:   
    Particle radius.

# Return
- pos::Matrix{Float64}  
    Particle positions.

- geometry_cfg::RectangleCfg  
    GeometryCfg that contains all particles.
"""
function rectangular_grid(num_p_x, num_p_y, offset, radius; NUM_T=Float64)
    x = Vector{NUM_T}()
    y = Vector{NUM_T}()
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
    
    geometry_cfg=RectangleCfg(
        length=num_p_x*2*radius + radius*offset*(num_p_x + 1),
        height=num_p_y*2*radius + radius*offset*(num_p_y + 1),
    )

    pos = Matrix{NUM_T}(undef, 2, num_p_x*num_p_y)
    pos[1, :] = x
    pos[2, :] = y
    return copy(reinterpret(SVector{size(pos, 1), NUM_T}, vec(pos))), geometry_cfg
end

"""
Create random velocitys for `num_p` particles with the maximum component 
value being `max_value`.
"""
function random_vel(num_p, max_value=1; rng=Random.default_rng(), NUM_T=Float64)
    vels = Vector{SVector{2, NUM_T}}(undef, num_p)
    for i in eachindex(vels)
        vels[i] = (rand(rng, SVector{2, NUM_T}) * 2 .- 1) * max_value
    end

    return vels
end

end