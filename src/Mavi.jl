module Mavi

include("configs.jl")

using .Configs

"Particles state (positions and velocities)"
@kwdef struct State{T}
    pos::Matrix{T}
    vel::Matrix{T}
end

include("space_checks.jl")
using .SpaceChecks

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
function System(;state::State{T}, space_cfg, dynamic_cfg, int_cfg) where {T}
    all_inside, out_ids = check_inside(state, space_cfg)
    if all_inside == false
        throw("Particles with ids=$(out_ids) outside space.")
    end
    num_p = size(state.pos)[2]
    diffs = Array{T, 3}(undef, 2, num_p, num_p)
    forces = Array{T, 2}(undef, 2, num_p)
    dists = zeros(T, num_p, num_p)
    System(state, space_cfg, dynamic_cfg, int_cfg, diffs, forces, dists, num_p)
end

include("integration.jl")
include("visualization.jl")
include("info.jl")

end
