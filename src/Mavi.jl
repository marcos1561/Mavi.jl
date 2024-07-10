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
Returns the indices of particles outside the given space configuration.
"""
function outside_particles(state::State, space_cfg::RectangleCfg)
    out_ids = findall((state.x .<= 0.0 .| state.x .>= space_cfg.length) .| (state.y .<= 0.0 .| state.y .>= space_cfg.height))
    return out_ids
end

function outside_particles(state::State, space_cfg::CircleCfg)
    r2 = state.x.^2 + state.y.^2
    out_ids = findall(r2 .>= space_cfg.radius^2)
    return out_ids
end
"""
Checks if all particles are inside the given space configuration. If not, throws an error.
"""
function check_inside(state::State, space_cfg::SpaceCfg)
    all_inside = true
    out_ids = outside_particles(state,space_cfg)
    if length(out_ids) != 0
        all_inside = false
    end
    return all_inside, out_ids
end

function System(;state::State{T}, space_cfg, dynamic_cfg, int_cfg) where {T}
    all_inside, out_ids = check_inside(state,space_cfg)
    if all_inside == false
        throw("particles $(out_ids) outside space.")
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
