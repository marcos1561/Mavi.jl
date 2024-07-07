module Mavi

include("configs.jl")

using .Configs

"Particles state (positions and velocities)"
@kwdef struct State{T}
    x::Vector{T}
    y::Vector{T}
    vx::Vector{T}
    vy::Vector{T}
end

"""
Base struct that represent a system of particles.

d=0: x axis
d=1: y axis
"""
struct System{T}
    state::State{T}
    space_cfg::SpaceCfg
    dynamic_cfg::DynamicCfg
    int_cfg::IntegrationCfg
    
    """
    Position difference between all particles
    
    diffs[d, i, j] = Position of particle i minus particles j in dimension d.

    d=1: x axis
    d=2: y axis
    """
    diffs::Array{T, 3}
    
    """
    Total force of all particles

    forces[d, i] = Total force on particle i in dimension d.
    """
    forces::Array{T, 2}
    
    """
    Distance between all particles

    dists[i, j] = Distance between particle i and j.
    """
    dists::Array{T, 2}

    "Number of particles"
    num_p::Int
end

function System(;state::State{T}, space_cfg, dynamic_cfg, int_cfg) where {T}
    num_p = length(state.x)
    
    diffs = Array{T, 3}(undef, 2, num_p, num_p)
    forces = Array{T, 2}(undef, 2, num_p)
    dists = zeros(T, num_p, num_p)
    System{T}(state, space_cfg, dynamic_cfg, int_cfg, diffs, forces, dists, num_p)
end

include("visualization.jl")

end
