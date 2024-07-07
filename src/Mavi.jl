module Mavi

include("configs.jl")

using .Configs

"Particles state (positions and velocities)"
@kwdef mutable struct State{T}
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
mutable struct System{T}
    state::State{T}
    space_cfg::SpaceCfg
    dynamic_cfg::DynamicCfg
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
    num_p = length(state.x)
    
    diffs = Array{T, 3}(undef, 2, num_p, num_p)
    forces = Array{T, 2}(undef, 2, num_p)
    dists = zeros(T, num_p, num_p)
    System{T}(state, space_cfg, dynamic_cfg, int_cfg, diffs, forces, dists, num_p)
end

"""
Calculates the total force acting on all particles

ko: coupling constant
ro: oscillation center
ra: maximum distance for which the potential is nonzero
"""
function forces!(system::System{T}) where {T}

    # Aliases
    state = system.state
    ko = system.dynamic_cfg.ko
    ro = system.dynamic_cfg.ro
    ra = system.dynamic_cfg.ra
    N = system.num_p
    # Initialize forces as zero
    system.forces .= 0.0
    for i in 1:N
        for j in i+1:N
            # Difference in position
            x_ij = state.x[i] - state.x[j]
            y_ij = state.y[i] - state.y[j]
            r_ij = sqrt(x_ij^2 + y_ij^2)
            # Check interaction range
            if r_ij < ra # cutoff distance
                fmod = -ko*(r_ij - ro) # restoring force
                fx_ij = fmod*x_ij/r_ij # unit vector x_ij/r_ij
                fy_ij = fmod*y_ij/r_ij
            else # too far
                fx_ij = 0.0
                fy_ij = 0.0
            end
            # Update values
            system.forces[1,i] += fx_ij
            system.forces[1,j] -= fx_ij
            system.forces[2,i] += fy_ij
            system.forces[2,j] -= fy_ij 
        end
    end
end

"""
Advances one time step using Velocity-Verlet

ko: coupling constant
ro: oscillation center
ra: maximum distance for which the potential is nonzero
"""
function step!(system::System{T}) where {T}

    # Aliases
    state = system.state
    N = system.num_p
    dt = system.int_cfg.dt # integration timestep
    # Truncated harmonic potential constants
    ko = system.dynamic_cfg.ko
    ro = system.dynamic_cfg.ro
    ra = system.dynamic_cfg.ra
    forces = system.forces
    diffs = system.diffs
    dists = system.dists

    # Calculate present forces
    forces!(system)
    old_forces = copy(system.forces)

    # Update positions
    term = dt^2/2 # quadratic term on accelerated movement
    state.x += state.vx*dt + forces[1,:]*term # uniformly accelerated movement
    state.y += state.vy*dt + forces[2,:]*term

    # Calculate new forces
    forces!(system)

    # Update velocities
    state.vx += dt*(system.forces[1,:] + old_forces[1,:])/2 # mean over old and new forces
    state.vy += dt*(system.forces[2,:] + old_forces[2,:])/2
end


end
