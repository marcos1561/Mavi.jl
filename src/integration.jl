module Integration

export newton_step!, szabo_step!, rtp_step!
export calc_forces!, calc_interaction, walls!
export calc_diff_and_dist, calc_diffs_and_dists!
export update_verlet!, update_rtp!, update_szabo!

using Base.Threads
using StaticArrays

using Mavi.Systems
using Mavi.States
using Mavi.Configs
using Mavi.ChunksMod

"Position difference (i - j) and distance between particle with id `i` and `j`"
@inline function calc_diff_and_dist(i, j, pos::Vector, space_cfg)
    @inbounds dr = pos[i] - pos[j]
    dist = sqrt(sum(abs2, dr))
    return dr, dist
end

@inline function calc_diff_and_dist(i, j, pos::Vector, space_cfg::SpaceCfg{PeriodicWalls, G}) where G <: RectangleCfg
    size_vec = space_cfg.geometry_cfg.size
    dr = pos[i] - pos[j]
    dr = dr - (abs.(dr) .> (size_vec / 2)) .* copysign.(size_vec, dr)
    dist = sqrt(sum(abs2, dr))
    return dr, dist
end

"Return force on particle with id `i` exerted by particle with id `j`."
function calc_interaction(i, j, dynamic_cfg::HarmTruncCfg, system::System)
    state = system.state
    space_cfg = system.space_cfg

    dr, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

    # Check interaction range
    if dist > dynamic_cfg.ra
        return SVector(0.0, 0.0)
    end

    fmod = -dynamic_cfg.ko*(dist - dynamic_cfg.ro) # restoring force
    # fx = fmod*x_ij/dist # unit vector x_ij/dist
    # fy = fmod*y_ij/dist

    return fmod/dist * dr
end

function calc_interaction(i, j, dynamic_cfg::LenJonesCfg, system::System)
    state = system.state
    space_cfg = system.space_cfg
    
    dr, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon

    # Force modulus
    fmod = 4*epsilon*(12*sigma^12/dist^13 - 6*sigma^6/dist^7)

    # Force components
    # fx = fmod*x_ij/dist
    # fy = fmod*y_ij/dist

    return fmod/dist * dr
end

function calc_interaction(i, j, dynamic_cfg::SzaboCfg, system::System)
    state = system.state
    space_cfg = system.space_cfg

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

function calc_interaction(i, j, dynamic_cfg::RunTumbleCfg, system::System)
    state = system.state
    space_cfg = system.space_cfg

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
function calc_forces!(system::System, chunks::Chunks, device::Sequencial)
    forces = get_forces(system)

    # Iteration over all chunks
    for col in 1:chunks.num_cols
        for row in 1:chunks.num_rows
            np = chunks.num_particles_in_chunk[row, col]
            # chunk = @view chunks.chunk_particles[1:np, row, col]
            chunk = @view chunks.chunk_particles[:, row, col]
            neighbors = chunks.neighbors[row, col]
            
            # Iteration over all particles in current chunk
            for i in 1:np
                p1_id = chunk[i]
                
                # Interaction between particles in the same chunk
                for j in (i+1):np
                    p2_id = chunk[j]
                    f1 = calc_interaction(p1_id, p2_id, system.dynamic_cfg, system)
                    
                    forces[p1_id] += f1
                    forces[p2_id] -= f1
                    # forces[1, p1_id] += fx
                    # forces[2, p1_id] += fy
                    # forces[1, p2_id] -= fx
                    # forces[2, p2_id] -= fy    
                end
                
                # Interaction of particles in neighboring chunks
                for neighbor_id in neighbors
                    nei_np = chunks.num_particles_in_chunk[neighbor_id]
                    # nei_chunk = @view chunks.chunk_particles[1:nei_np, neighbor_id]
                    nei_chunk = @view chunks.chunk_particles[:, neighbor_id]
                    
                    # println("num: ", nei_np)
                    # println("vec: ", chunks.chunk_particles[:, neighbor_id])
                    # println("vec slice: ", chunks.chunk_particles[1:nei_np, neighbor_id])
                    # println("====")

                    for j in 1:nei_np
                        p2_id = nei_chunk[j]
                        f1 = calc_interaction(p1_id, nei_chunk[j], 
                            system.dynamic_cfg, system)
                        
                        forces[p1_id] +=  f1
                        forces[p2_id] -=  f1
                        # forces[1, p1_id] +=  fx
                        # forces[2, p1_id] +=  fy
                        # forces[1, p2_id] -= fx
                        # forces[2, p2_id] -= fy        
                    end
                end
            end
        end
    end
end

function calc_forces!(system::System, chunks::Chunks, device::Threaded)
    # forces_local = [zeros(size(system.forces)) for _ in 1:nthreads()]

    Threads.@threads :static for col in 1:chunks.num_cols
        # forces = @view system.forces[:, :, Threads.threadid()]
        forces = system.forces[Threads.threadid()]
        for row in 1:chunks.num_rows
            np = chunks.num_particles_in_chunk[row, col]
            chunk = @view chunks.chunk_particles[:, row, col]
            neighbors = chunks.neighbors[row, col]
            
            for i in 1:np
                p1_id = chunk[i]
                for j in (i+1):np
                    p2_id = chunk[j]
                    f1 = calc_interaction(p1_id, p2_id, 
                        system.dynamic_cfg, system)
                    
                    forces[p1_id] += f1
                    forces[p2_id] -= f1
                    # forces[1, p1_id] += fx
                    # forces[2, p1_id] += fy
                    # forces[1, p2_id] -= fx
                    # forces[2, p2_id] -= fy    
                end
                for neighbor_id in neighbors
                    nei_np = chunks.num_particles_in_chunk[neighbor_id]
                    nei_chunk = @view chunks.chunk_particles[:, neighbor_id]
                    for j in 1:nei_np
                        p2_id = nei_chunk[j]
                        f1 = calc_interaction(p1_id, nei_chunk[j], 
                            system.dynamic_cfg, system)
                        forces[p1_id] +=  f1
                        forces[p2_id] -=  f1
                        # forces[1, p1_id] +=  fx
                        # forces[2, p1_id] +=  fy
                        # forces[1, p2_id] -= fx
                        # forces[2, p2_id] -= fy        
                    end
                end
            end
        end
    end

    # Soma todas as forças locais no array global
    # get_forces(system) .= sum(system.forces, dims=3)
    get_forces(system) .= sum(system.forces)
end
# function calc_forces!(system::System, chunks::Chunks, device::Threaded)
#     # chunk_locks = Vector{Threads.ReentrantLock}(undef, chunks.num_cols * chunks.num_rows)
#     chunk_locks = Matrix{Threads.ReentrantLock}(undef, chunks.num_rows, chunks.num_cols)
#     for i in 1:chunks.num_rows
#         for j in 1:chunks.num_cols
#             chunk_locks[i, j] = Threads.ReentrantLock()
#         end
#     end

#     forces = system.forces 
#     forces .= 0

#     # @threads for col in 1:chunks.num_cols
#     #     for row in 1:chunks.num_rows
#     @threads for cart_idx in CartesianIndices((chunks.num_rows, chunks.num_cols))
#             np = chunks.num_particles_in_chunk[cart_idx]
#             chunk = @view chunks.chunk_particles[1:np, cart_idx]
#             neighbors = chunks.neighbors[cart_idx]
            
#             for i in 1:np
#                 p1_id = chunk[i]
#                 for j in (i+1):np
#                     p2_id = chunk[j]
#                     fx, fy = calc_interaction(p1_id, p2_id, 
#                         system.state, system.dynamic_cfg, system.space_cfg)
                    
#                     lock(chunk_locks[cart_idx]) do
#                         forces[1, p1_id] += fx
#                         forces[2, p1_id] += fy
#                         forces[1, p2_id] -= fx
#                         forces[2, p2_id] -= fy    
#                     end
#                 end
#                 for neighbor_id in neighbors
#                     nei_np = chunks.num_particles_in_chunk[neighbor_id]
#                     nei_chunk = @view chunks.chunk_particles[1:nei_np, neighbor_id]
#                     for j in 1:nei_np
#                         p2_id = nei_chunk[j]
#                         fx, fy = calc_interaction(p1_id, nei_chunk[j], 
#                             system.state, system.dynamic_cfg, system.space_cfg)
                        
#                         lock(chunk_locks[cart_idx]) do
#                             forces[1, p1_id] +=  fx
#                             forces[2, p1_id] +=  fy
#                         end
#                         lock(chunk_locks[neighbor_id]) do
#                             forces[1, p2_id] -= fx
#                             forces[2, p2_id] -= fy
#                         end        
#                     end
#                 end
#             end
#         end
#     # end
# end

"Compute total forces on particles."
function calc_forces!(system::System, chunks::Nothing, device::Sequencial)
    forces = get_forces(system)

    N = system.num_p
    for i in 1:N
        for j in i+1:N
            if !is_valid_pair(system.state, system.dynamic_cfg, i, j)
                continue
            end
            f = calc_interaction(i, j, system.dynamic_cfg, system)
            forces[i] += f
            forces[j] -= f
        end
    end
end

@inline calc_forces!(system::System) = calc_forces!(system, system.chunks, system.int_cfg.device)

"Rigid walls collisions. Reflect velocity on collision."
function walls!(system::System, space_cfg::SpaceCfg{RigidWalls, G}) where G <: RectangleCfg
    state = system.state
    geometry_cfg = space_cfg.geometry_cfg
    size_vec = SVector(geometry_cfg.length, geometry_cfg.height)
    r = particle_radius(system.dynamic_cfg)

    for i in eachindex(state.pos)
        rel_pos = state.pos[i] - geometry_cfg.bottom_left

        out = (rel_pos .+ r) .> size_vec .|| (rel_pos .- r) .< 0
        if any(out)
            state.vel[i] = state.vel[i] .* (-2*out .+ 1)
        end

        # if rel_pos[1] + r > geometry_cfg.length || rel_pos[1] - r < 0 
        #     state.vel[i][1] .*= -1
        # end
        # if rel_pos[2] + r > geometry_cfg.height || rel_pos[2] - r < 0 
        #     state.vel[i][2] .*= -1
        # end
    end

    # x = state.pos[1, :] .- geometry_cfg.bottom_left[1]
    # y = state.pos[2, :] .- geometry_cfg.bottom_left[2]

    # out_x = @. ((x + r) > geometry_cfg.length) || ((x - r) < 0)
    # out_y = @. ((y + r) > geometry_cfg.height) || ((y - r) < 0)
    
    # # out_x = @. ((state.pos[1, :] - geometry_cfg.bottom_left[1] + r) > geometry_cfg.length) || ((state.pos[1, :] - geometry_cfg.bottom_left[1] - r) < 0)
    # # out_y = @. ((state.pos[2, :] - geometry_cfg.bottom_left[2] + r) > geometry_cfg.height) || ((state.pos[2, :] - geometry_cfg.bottom_left[2] - r) < 0)

    # state.vel[1, findall(out_x)] .*= -1.0
    # state.vel[2, findall(out_y)] .*= -1.0
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

"Periodic walls"
function walls!(system::System, space_cfg::SpaceCfg{PeriodicWalls, G}) where G <: RectangleCfg
    state = system.state
    geometry_cfg = space_cfg.geometry_cfg
    # center = (geometry_cfg.bottom_left[1] + geometry_cfg.length/2, geometry_cfg.bottom_left[2] + geometry_cfg.height/2)
    half_size = SVector(geometry_cfg.length, geometry_cfg.height) / 2
    center = geometry_cfg.bottom_left + half_size 
    
    for i in eachindex(state.pos)
        pos = state.pos[i]
        diff = pos - center
        
        out_bounds = abs.(diff) .> half_size
        if any(out_bounds)
            state.pos[i] = pos .- sign.(diff) .* (half_size * 2) .* out_bounds
        end
    end

    # context = (center, geometry_cfg)
    # iterate_particles(state, context, resolve_wall!)

    # for i in 1:system.num_p
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
end

@inline walls!(system::System) = walls!(system, system.space_cfg, system.dynamic_cfg)
@inline walls!(system::System, space_cfg::SpaceCfg, dynamic_cfg::DynamicCfg) = walls!(system, system.space_cfg)

"Update state using Velocity-Verlet"
function update_verlet!(system::System, calc_forces_in!)
    state = system.state
    forces = get_forces(system)
    # old_forces = copy(forces)
    old_forces = forces

    dt = system.int_cfg.dt

    # Update positions
    term = dt^2/2 # quadratic term on accelerated movement
    @. state.pos += state.vel * dt + forces * term

    clean_forces!(system)
    calc_forces_in!(system)

    # Update velocities
    @. state.vel += dt/2 * (forces + old_forces)
end

function update_szabo!(system::System)
    state::SelfPropelledState = system.state
    forces = get_forces(system)
    dt = system.int_cfg.dt

    dynamic_cfg = system.dynamic_cfg
    vo, relax_time = dynamic_cfg.vo, dynamic_cfg.relax_time
    mu, dr = dynamic_cfg.mobility, dynamic_cfg.rot_diff 

    for i in 1:system.num_p
        theta = state.pol_angle[i]
        pol_x, pol_y = cos(theta), sin(theta) 

        vel_x = vo * pol_x + mu * forces[1, i]
        vel_y = vo * pol_y + mu * forces[2, i]
        
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
    forces = get_forces(system)

    dt = system.int_cfg.dt

    dynamic_cfg = system.dynamic_cfg
    vo, sigma, epsilon, tumble_rate = dynamic_cfg.vo, dynamic_cfg.sigma, dynamic_cfg.epsilon, dynamic_cfg.tumble_rate

    for i in 1:system.num_p
        theta = state.pol_angle[i]
        pol_x, pol_y = cos(theta), sin(theta) 

        vel_x = vo * pol_x + forces[1, i]
        vel_y = vo * pol_y + forces[2, i]
        
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
function newton_step!(system::System)
    clean_forces!(system)
    update_chunks!(system.chunks)
    calc_forces!(system)
    update_verlet!(system, calc_forces!)
    walls!(system, system.space_cfg)
end

function szabo_step!(system::System)
    clean_forces!(system)
    update_chunks!(system.chunks)
    calc_forces!(system)
    update_szabo!(system)
    walls!(system, system.space_cfg)
end

function rtp_step!(system::System)
    clean_forces!(system)
    update_chunks!(system.chunks)
    calc_forces!(system)
    update_rtp!(system)
    walls!(system, system.space_cfg)
end

end