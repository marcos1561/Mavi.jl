module Integration

using Mavi
using Mavi.Configs: IntCfg, ChunksIntCfg
using Mavi.Configs
using Mavi.ChunksMod

"Position difference and distance between particle with id `i` and `j`"
function calc_diff_and_dist(i, j, state::State)
    dx = state.pos[1, i] - state.pos[1, j]
    dy = state.pos[2, i] - state.pos[2, j]
    dist = sqrt(dx^2 + dy^2)
    return dx, dy, dist
end

"Compute the difference and distance between all particles."
function calc_diffs_and_dists!(system::System)
    state = system.state

    for i in 1:system.num_p
        for j in i+1:system.num_p
            x_ij, y_ij, dist = calc_diff_and_dist(i, j, state)

            # Update values
            system.diffs[1,i,j] = x_ij
            system.diffs[1,j,i] = -x_ij
            system.diffs[2,i,j] = y_ij
            system.diffs[2,j,i] = -y_ij
            system.dists[i,j] = dist
            system.dists[j,i] = dist
        end
    end
end

"Return force on particle with id `i` exerted by particle with id `j`."
function calc_interaction(i, j, state::State, dynamic_cfg::HarmTruncCfg)
    x_ij, y_ij, dist = calc_diff_and_dist(i, j, state)

    # Check interaction range
    if dist > dynamic_cfg.ra
        return 0.0, 0.0
    end

    fmod = -dynamic_cfg.ko*(dist - dynamic_cfg.ro) # restoring force
    fx = fmod*x_ij/dist # unit vector x_ij/dist
    fy = fmod*y_ij/dist

    return fx, fy
end

function calc_interaction(i, j, state::State, dynamic_cfg::LenJonesCfg)
    x_ij, y_ij, dist = calc_diff_and_dist(i, j, state)

    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon

    # Force modulus
    fmod = 4*epsilon*(12*sigma^12/dist^13 - 6*sigma^6/dist^7)

    # Force components
    fx = fmod*x_ij/dist
    fy = fmod*y_ij/dist

    return fx, fy
end

"Compute total forces on particles using chucks."
function calc_forces_chunks!(system::System)
    system.forces .= 0
    
    chunks = system.chunks

    # Iteration over all chunks
    for col in 1:chunks.num_cols
        for row in 1:chunks.num_rows
            np = chunks.num_particles_in_chunk[row, col]
            chunk = @view chunks.chunk_particles[1:np, row, col]
            neighbors = chunks.neighbors[row, col]
            
            # Iteration over all particles in current chunk
            for i in 1:np
                p1_id = chunk[i]
                
                # Interaction between particles in the same chunk
                for j in (i+1):np
                    p2_id = chunk[j]
                    fx, fy = calc_interaction(p1_id, p2_id, 
                        system.state, system.dynamic_cfg)
                    
                    system.forces[1, p1_id] += fx
                    system.forces[2, p1_id] += fy
                    system.forces[1, p2_id] -= fx
                    system.forces[2, p2_id] -= fy    
                end
                
                # Interaction of particles in neighboring chunks
                for neighbor_id in neighbors
                    nei_np = chunks.num_particles_in_chunk[neighbor_id]
                    nei_chunk = @view chunks.chunk_particles[1:nei_np, neighbor_id]
                    
                    for j in 1:nei_np
                        p2_id = nei_chunk[j]
                        fx, fy = calc_interaction(p1_id, nei_chunk[j], 
                            system.state, system.dynamic_cfg)
                        
                        system.forces[1, p1_id] +=  fx
                        system.forces[2, p1_id] +=  fy
                        system.forces[1, p2_id] -= fx
                        system.forces[2, p2_id] -= fy        
                    end
                end
            end
        end
    end
end

"Compute total forces on particles."
function calc_forces!(system::System)
    system.forces .= 0

    N = system.num_p
    for i in 1:N
        for j in i+1:N
            fx, fy = calc_interaction(i, j, system.state, system.dynamic_cfg)
            
            system.forces[1, i] +=  fx
            system.forces[2, i] +=  fy
            system.forces[1, j] -= fx
            system.forces[2, j] -= fy   
        end
    end
end

"Rigid walls collisions. Reflect velocity on collision."
function rigid_walls!(system::System, space_cfg::RectangleCfg)
    state = system.state
    r = particle_radius(system.dynamic_cfg)

    x = state.pos[1, :] .- space_cfg.bottom_left[1]
    y = state.pos[2, :] .- space_cfg.bottom_left[2]

    out_x = @. ((x + r) > space_cfg.length) || ((x - r) < 0)
    out_y = @. ((y + r) > space_cfg.height) || ((y - r) < 0)
    
    # out_x = @. ((state.pos[1, :] - space_cfg.bottom_left[1] + r) > space_cfg.length) || ((state.pos[1, :] - space_cfg.bottom_left[1] - r) < 0)
    # out_y = @. ((state.pos[2, :] - space_cfg.bottom_left[2] + r) > space_cfg.height) || ((state.pos[2, :] - space_cfg.bottom_left[2] - r) < 0)

    state.vel[1, findall(out_x)] .*= -1.0
    state.vel[2, findall(out_y)] .*= -1.0
end

function rigid_walls!(system::System, space_cfg::CircleCfg)
    state = system.state
    max_r2 = (space_cfg.radius - particle_radius(system.dynamic_cfg))^2
    for i in 1:system.num_p
        x, y = state.pos[1, i], state.pos[2, i]
        x2, y2 = x^2, y^2
        if (x2 + y2) > max_r2
            r2 = x2 + y2
            vx, vy = state.vel[1, i], state.vel[2, i]
            state.vel[1, i] = (vx*(y2 - x2) - 2*vy*x*y) / r2
            state.vel[2, i] = (-vy*(y2 - x2) - 2*vx*x*y) / r2
        end
    end
end

"Update state using Velocity-Verlet"
function update_verlet!(system::System, calc_forces_in!)
    state = system.state
    old_forces = copy(system.forces)

    dt = system.int_cfg.dt

    # Update positions
    term = dt^2/2 # quadratic term on accelerated movement
    state.pos .+= state.vel * dt + system.forces * term

    calc_forces_in!(system)

    # Update velocities
    state.vel .+= dt/2 * (system.forces + old_forces)
end

"Advance system one time step."
function step!(system::System, int_cfg::IntCfg)
    calc_forces!(system)
    update_verlet!(system, calc_forces!)
    rigid_walls!(system, system.space_cfg)
end

function step!(system::System, int_cfg::ChunksIntCfg)
    update_chunks!(system.chunks)
    calc_forces_chunks!(system)
    update_verlet!(system, calc_forces_chunks!)
    rigid_walls!(system, system.space_cfg)
end

end