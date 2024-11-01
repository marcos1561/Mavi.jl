module Utils

export get_ring_radius

function get_ring_radius(p_radius, num_particles)
    return (p_radius*2) / (2 * (1 - cos(2*π/num_particles)))^.5
end

end