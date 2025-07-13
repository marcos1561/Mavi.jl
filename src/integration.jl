module Integration

using Mavi.Systems: System
using Mavi.States
using Mavi.Configs
using Mavi.ChunksMod

"Position difference (i - j) and distance between particle with id `i` and `j`"
function calc_diff_and_dist(i, j, pos, space_cfg)
    dx = pos[1, i] - pos[1, j]
    dy = pos[2, i] - pos[2, j]
    dist = sqrt(dx^2 + dy^2)
    return dx, dy, dist
end

function calc_diff_and_dist(i, j, pos, space_cfg::SpaceCfg{PeriodicWalls, RectangleCfg})
    dx = pos[1, i] - pos[1, j]
    dy = pos[2, i] - pos[2, j]
    
    geometry_cfg = space_cfg.geometry_cfg
    length, height = geometry_cfg.length, geometry_cfg.height 
    
    if (abs(dx) > length * 0.5)
        dx -= copysign(length, dx)
    end

    if (abs(dy) > height * 0.5)
        dy -= copysign(height, dy)
    end

    dist = sqrt(dx^2 + dy^2)
    return dx, dy, dist
end

"Compute the difference and distance between all particles."
function calc_diffs_and_dists!(system::System, space_cfg)
    state = system.state

    for i in 1:system.num_p
        for j in i+1:system.num_p
            x_ij, y_ij, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

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
function calc_interaction(i, j, state::State, dynamic_cfg::HarmTruncCfg, space_cfg)
    x_ij, y_ij, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

    # Check interaction range
    if dist > dynamic_cfg.ra
        return 0.0, 0.0
    end

    fmod = -dynamic_cfg.ko*(dist - dynamic_cfg.ro) # restoring force
    fx = fmod*x_ij/dist # unit vector x_ij/dist
    fy = fmod*y_ij/dist

    return fx, fy
end

function calc_interaction(i, j, state::State, dynamic_cfg::LenJonesCfg, space_cfg)
    x_ij, y_ij, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon

    # Force modulus
    fmod = 4*epsilon*(12*sigma^12/dist^13 - 6*sigma^6/dist^7)

    # Force components
    fx = fmod*x_ij/dist
    fy = fmod*y_ij/dist

    return fx, fy
end

function calc_interaction(i, j, state::State, dynamic_cfg::SzaboCfg, space_cfg)
    dx, dy, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

    if dist > dynamic_cfg.r_max
        return 0, 0
    end

    if dist > dynamic_cfg.r_eq
        f_mod = dynamic_cfg.k_adh / (dynamic_cfg.r_eq)
    else
        f_mod = dynamic_cfg.k_rep / (dynamic_cfg.r_max - dynamic_cfg.r_eq)
    end

    r = dist - dynamic_cfg.r_eq
    fx, fy = -f_mod * dx * r, -f_mod * dy * r
    return fx, fy
end

function calc_interaction(i, j, state::State, dynamic_cfg::RunTumbleCfg, space_cfg)
    x_ij, y_ij, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon
    cutoff = 2^(1/6)*sigma

    # Distant particles
    if dist > cutoff
        return 0.0, 0.0
    end

    # Force modulus
    fmod = -4*epsilon*(-12*sigma^12/dist^13 + 6*sigma^6/dist^7)

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
                        system.state, system.dynamic_cfg, system.space_cfg)
                    
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
                            system.state, system.dynamic_cfg, system.space_cfg)
                        
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
            fx, fy = calc_interaction(i, j, system.state, system.dynamic_cfg, system.space_cfg)
            
            system.forces[1, i] +=  fx
            system.forces[2, i] +=  fy
            system.forces[1, j] -= fx
            system.forces[2, j] -= fy   
        end
    end
end

"Rigid walls collisions. Reflect velocity on collision."
function walls!(system::System, space_cfg::SpaceCfg{RigidWalls, RectangleCfg})
    state = system.state
    geometry_cfg = space_cfg.geometry_cfg
    r = particle_radius(system.dynamic_cfg)

    x = state.pos[1, :] .- geometry_cfg.bottom_left[1]
    y = state.pos[2, :] .- geometry_cfg.bottom_left[2]

    out_x = @. ((x + r) > geometry_cfg.length) || ((x - r) < 0)
    out_y = @. ((y + r) > geometry_cfg.height) || ((y - r) < 0)
    
    # out_x = @. ((state.pos[1, :] - geometry_cfg.bottom_left[1] + r) > geometry_cfg.length) || ((state.pos[1, :] - geometry_cfg.bottom_left[1] - r) < 0)
    # out_y = @. ((state.pos[2, :] - geometry_cfg.bottom_left[2] + r) > geometry_cfg.height) || ((state.pos[2, :] - geometry_cfg.bottom_left[2] - r) < 0)

    state.vel[1, findall(out_x)] .*= -1.0
    state.vel[2, findall(out_y)] .*= -1.0
end

function walls!(system::System, space_cfg::SpaceCfg{RigidWalls, CircleCfg})
    state = system.state
    max_r2 = (space_cfg.geometry_cfg.radius - particle_radius(system.dynamic_cfg))^2
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

# function resolve_wall!(center, geometry_cfg, state)
#     pos_i = @view state.pos[:, i] 
#     diff = pos_i[1] - center[1]  
#     if abs(diff) > geometry_cfg.length/2
#         pos_i[1] -= sign(diff) * geometry_cfg.length
#     end

#     diff = pos_i[2] - center[2]  
#     if abs(diff) > geometry_cfg.height/2
#         pos_i[2] -= sign(diff) * geometry_cfg.height
#     end
# end

"Periodic walls"
function walls!(system::System, space_cfg::SpaceCfg{PeriodicWalls, RectangleCfg}) 
    state = system.state
    geometry_cfg = space_cfg.geometry_cfg
    center = (geometry_cfg.bottom_left[1] + geometry_cfg.length/2, geometry_cfg.bottom_left[2] + geometry_cfg.height/2)

    # context = (center, geometry_cfg)
    # iterate_particles(state, context, resolve_wall!)

    for i in 1:system.num_p
        pos_i = @view state.pos[:, i] 
        diff = pos_i[1] - center[1]  
        if abs(diff) > geometry_cfg.length/2
            pos_i[1] -= sign(diff) * geometry_cfg.length
        end

        diff = pos_i[2] - center[2]  
        if abs(diff) > geometry_cfg.height/2
            pos_i[2] -= sign(diff) * geometry_cfg.height
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

function update_szabo!(system::System)
    state::SelfPropelledState = system.state
    
    dt = system.int_cfg.dt

    dynamic_cfg = system.dynamic_cfg
    vo, relax_time = dynamic_cfg.vo, dynamic_cfg.relax_time
    mu, dr = dynamic_cfg.mobility, dynamic_cfg.rot_diff 

    for i in 1:system.num_p
        theta = state.pol_angle[i]
        pol_x, pol_y = cos(theta), sin(theta) 

        vel_x = vo * pol_x + mu * system.forces[1, i]
        vel_y = vo * pol_y + mu * system.forces[2, i]
        
        speed = sqrt(vel_x^2 + vel_y^2)
    
        cross_prod = (pol_x * vel_y - pol_y * vel_x) / speed 
        if abs(cross_prod) > 1
            cross_prod = sign(cross_prod)
        end

        d_theta = 1/relax_time * asin(cross_prod) * dt + sqrt(2 * dr * dt) * randn()
        
        state.pos[1, i] += vel_x * dt
        state.pos[2, i] += vel_y * dt
        state.pol_angle[i] += d_theta
    end
end

function update_rtp!(system::System)
    state::SelfPropelledState = system.state
    
    dt = system.int_cfg.dt

    dynamic_cfg = system.dynamic_cfg
    vo, sigma, epsilon, tumble_rate = dynamic_cfg.vo, dynamic_cfg.sigma, dynamic_cfg.epsilon, dynamic_cfg.tumble_rate

    for i in 1:system.num_p
        theta = state.pol_angle[i]
        pol_x, pol_y = cos(theta), sin(theta) 

        vel_x = vo * pol_x + system.forces[1, i]
        vel_y = vo * pol_y + system.forces[2, i]
        
        speed = sqrt(vel_x^2 + vel_y^2)
        
        # Update positions
        state.pos[1, i] += vel_x * dt
        state.pos[2, i] += vel_y * dt

        # Update director
        u = rand()
        if u < tumble_rate * dt # accept tumble
            state.pol_angle[i] = 2*π*rand()
        end
    end
end

"Advance system one time step."
function step!(system::System, int_cfg::IntCfg)
    calc_forces!(system)
    update_verlet!(system, calc_forces!)
    walls!(system, system.space_cfg)
end

function step!(system::System, int_cfg::ChunksIntCfg)
    update_chunks!(system.chunks)
    calc_forces_chunks!(system)
    update_verlet!(system, calc_forces_chunks!)
    walls!(system, system.space_cfg)
end

function szabo_step!(system::System, int_cfg::IntCfg)
    calc_forces!(system)
    update_szabo!(system)
    walls!(system, system.space_cfg)
end

function szabo_step!(system::System, int_cfg::ChunksIntCfg)
    update_chunks!(system.chunks)
    calc_forces_chunks!(system)
    update_szabo!(system)
    walls!(system, system.space_cfg)
end

function rtp_step!(system::System, int_cfg::IntCfg)
    calc_forces!(system)
    update_rtp!(system)
    walls!(system, system.space_cfg)
end

function rtp_step!(system::System, int_cfg::ChunksIntCfg)
    update_chunks!(system.chunks)
    calc_forces_chunks!(system)
    update_rtp!(system)
    walls!(system, system.space_cfg)
end

end