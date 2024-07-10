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

d=1: x axis
d=2: y axis
"""
struct System{T, C1<:SpaceCfg, C2<:DynamicCfg}
    state::State{T}
    space_cfg::C1
    dynamic_cfg::C2
    int_cfg::IntegrationCfg
    
    """
    Position difference between all particles
    
    diffs[d, i, j] = Position of particle i minus particle j in dimension d.

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

"""
Checks if all particles are inside the given configuration space. If not, throws an error.
"""
function check_inside(state::State, space_cfg::RectangleCfg)
    if all(x -> 0.0 <= x <= space_cfg.length, state.x) && all(y -> 0.0 <= y <= space_cfg.height, state.y)
        return true
    else
        return false
    end
end
function check_inside(state::State{T}, space_cfg::CircleCfg) where {T}
    r2 = state.x.^2 + state.y.^2
    if all(r -> r <= space_cfg.radius^2, r2)
        return true
    else
        return false
    end
end


function System(;state::State{T}, space_cfg, dynamic_cfg, int_cfg) where {T}
    if check_inside(state,space_cfg) == false
        throw("particles must be inside space")
    end
    num_p = length(state.x)
    diffs = Array{T, 3}(undef, 2, num_p, num_p)
    forces = Array{T, 2}(undef, 2, num_p)
    dists = zeros(T, num_p, num_p)
    System(state, space_cfg, dynamic_cfg, int_cfg, diffs, forces, dists, num_p)
end

include("integration.jl")
include("visualization.jl")

end
