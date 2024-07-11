module Integration

using Mavi: System
using Mavi.Configs

"Compute the difference and distance between all particles."
function calc_diffs_and_dists!(system::System)
    state = system.state

    for i in 1:system.num_p
        for j in i+1:system.num_p
            # Difference in position
            x_ij = state.pos[1, i] - state.pos[1, j]
            y_ij = state.pos[2, i] - state.pos[2, j]
            
            r_ij = sqrt(x_ij^2 + y_ij^2)

            # Update values
            system.diffs[1,i,j] = x_ij
            system.diffs[1,j,i] = -x_ij
            system.diffs[2,i,j] = y_ij
            system.diffs[2,j,i] = -y_ij
            system.dists[i,j] = r_ij
            system.dists[j,i] = r_ij
        end
    end
end

"""
Compute total force acting on all particles using
a truncated harmonic potencial.
"""
function calc_forces!(system::System, dynamic_cfg::HarmTruncCfg)
    # Aliases
    dists = system.dists
    diffs = system.diffs

    ko = dynamic_cfg.ko
    ro = dynamic_cfg.ro
    ra = dynamic_cfg.ra
    N = system.num_p
    # Initialize forces as zero
    system.forces .= 0.0
    for i in 1:N
        for j in i+1:N
            dist = dists[i, j]

            # Check interaction range
            if dist < ra # cutoff distance
                fmod = -ko*(dist - ro) # restoring force
                fx_ij = fmod*diffs[1, i, j]/dist # unit vector x_ij/dist
                fy_ij = fmod*diffs[2, i, j]/dist
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

"Rigid walls collisions. Reflect velocity on collision."
function rigid_walls!(system::System, space_cfg::RectangleCfg)
    state = system.state
    r = particle_radius(system.dynamic_cfg)

    out_x = @. ((state.pos[1, :] + r) > space_cfg.length) || ((state.pos[1, :] - r) < 0)
    out_y = @. ((state.pos[2, :] + r) > space_cfg.height) || ((state.pos[2, :] - r) < 0)

    state.vel[1, findall(out_x)] .*= -1.0
    state.vel[2, findall(out_y)] .*= -1.0

    # for i in 1:system.num_p
    #     # Hit vertical wall
    #     if ((state.x[i]+r) > space_cfg.length) || ((state.x[i]-r) < 0)
    #         state.vx[i] *= -1.
    #     end

    #     # Hit horizontal wall
    #     if ((state.y[i]+r) > space_cfg.height) || ((state.y[i]-r) < 0)
    #         state.vy[i] *= -1.
    #     end
    # end
end

function rigid_walls!(system::System, space_cfg::CircleCfg)
    state = system.state
    max_r2 = (space_cfg.radius - particle_radius(system.dynamic_cfg))^2
    for i in 1:system.num_p
        x, y = state.x[i], state.y[i]
        x2, y2 = x^2, y^2
        if (x2 + y2) > max_r2
            r2 = x2 + y2
            vx, vy = state.vx[i], state.vy[i]
            state.vx[i] = (vx*(y2 - x2) - 2*vy*x*y) / r2
            state.vy[i] = (-vy*(y2 - x2) - 2*vx*x*y) / r2
        end
    end
end

"Update state using Velocity-Verlet"
function update_verlet!(system::System)
    state = system.state
    old_forces = copy(system.forces)

    dt = system.int_cfg.dt

    # Update positions
    term = dt^2/2 # quadratic term on accelerated movement
    # state.x .+= state.vx*dt + system.forces[1,:]*term # uniformly accelerated movement
    # state.y .+= state.vy*dt + system.forces[2,:]*term
    state.pos .+= state.vel * dt + system.forces * term

    # Calculate new forces
    calc_diffs_and_dists!(system)
    calc_forces!(system, system.dynamic_cfg)

    # Update velocities
    # state.vx .+= dt*(system.forces[1,:] + old_forces[1,:])/2 # mean over old and new forces
    # state.vy .+= dt*(system.forces[2,:] + old_forces[2,:])/2
    state.vel .+= dt/2 * (system.forces + old_forces)
end

"Advance system one time step."
function step!(system::System)
    calc_diffs_and_dists!(system)
    calc_forces!(system, system.dynamic_cfg)
    update_verlet!(system)
    rigid_walls!(system, system.space_cfg)
end

end