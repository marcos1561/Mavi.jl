module Integration

export newton_step!, szabo_step!, rtp_step!
export calc_forces!, calc_interaction, walls!
export calc_diff, calc_diffs_and_dists!
export update_verlet!, update_rtp!, update_szabo!
export get_step_function

using Base.Threads
using StaticArrays

using Mavi.Systems
using Mavi.States
using Mavi.Configs
using Mavi.ChunksMod

"Position difference (r1 - r2)"
# @inline function calc_diff(i, j, pos::Vector, space_cfg)
@inline function calc_diff(r1, r2, space_cfg)
    @inbounds dr = r1 - r2
    return dr
end

@inline function calc_diff(r1, r2, space_cfg::SpaceCfg{PeriodicWalls, G}) where G <: RectangleCfg
    size_vec = space_cfg.geometry_cfg.size
    dr = r1 - r2
    dr = dr - (abs.(dr) .> (size_vec / 2)) .* copysign.(size_vec, dr)
    return dr
end

@inline function calc_diff(r1, r2, space_cfg::SpaceCfg{W, G}) where {W <: ManyWalls, G <: ManyGeometries}
    calc_diff(r1, r2, SpaceCfg(space_cfg.wall_type.list[1], space_cfg.geometry_cfg.list[1]))
end

"Return force on particle with id `i` exerted by particle with id `j`."
function calc_interaction(i, j, dynamic_cfg::HarmTruncCfg, system::System)
    pos = system.state.pos
    space_cfg = system.space_cfg

    dr = calc_diff(pos[i], pos[j], space_cfg)
    dist = sqrt(sum(dr.^2))

    # Check interaction range
    if dist > dynamic_cfg.ra
        return zero(dr)
    end

    fmod = -dynamic_cfg.ko*(dist - dynamic_cfg.ro) # restoring force

    return fmod/dist * dr
end

function calc_interaction(i, j, dynamic_cfg::LenJonesCfg, system::System)
    pos = system.state.pos
    space_cfg = system.space_cfg
    
    dr = calc_diff(pos[i], pos[j], space_cfg)
    dist = sqrt(sum(dr.^2))

    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon

    # Force modulus
    fmod = 4*epsilon*(12*sigma^12/dist^13 - 6*sigma^6/dist^7)

    return fmod/dist * dr
end

function calc_interaction(i, j, dynamic_cfg::SzaboCfg, system::System)
    pos = system.state.pos
    space_cfg = system.space_cfg

    dr = calc_diff(pos[i], pos[j], space_cfg)
    dist = sqrt(sum(dr.^2))

    if dist > dynamic_cfg.r_max
        return zero(dr)
    end

    if dist > dynamic_cfg.r_eq
        f_mod = dynamic_cfg.k_adh / (dynamic_cfg.r_eq)
    else
        f_mod = dynamic_cfg.k_rep / (dynamic_cfg.r_max - dynamic_cfg.r_eq)
    end

    r = dist - dynamic_cfg.r_eq
    return -f_mod * r * dr
end

function calc_interaction(i, j, dynamic_cfg::RunTumbleCfg, system::System)
    pos = system.state.pos
    space_cfg = system.space_cfg

    dr = calc_diff(pos[i], pos[j], space_cfg)
    dist = sqrt(sum(dr.^2))

    sigma = dynamic_cfg.sigma
    epsilon = dynamic_cfg.epsilon
    cutoff = 2^(1/6)*sigma

    # Distant particles
    if dist > cutoff
        return zero(dr)
    end

    # Force modulus
    fmod = -4*epsilon*(-12*sigma^12/dist^13 + 6*sigma^6/dist^7)

    return fmod / dist * dr 
end

"Compute total forces on particles using chucks."
function calc_forces!(system::System, chunks::Chunks, device::Sequencial)
    forces = get_forces(system)

    # Iteration over all chunks
    for col in 1:chunks.num_cols
        for row in 1:chunks.num_rows
            np = chunks.num_particles_in_chunk[row, col]
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
                end
                
                # Interaction of particles in neighboring chunks
                for neighbor_id in neighbors
                    nei_np = chunks.num_particles_in_chunk[neighbor_id]
                    nei_chunk = @view chunks.chunk_particles[:, neighbor_id]
                    
                    # println("num: ", nei_np)
                    # println("vec: ", chunks.chunk_particles[:, neighbor_id])
                    # println("vec slice: ", chunks.chunk_particles[1:nei_np, neighbor_id])
                    # println("====")

                    for j in 1:nei_np
                        p2_id = nei_chunk[j]
                        f1 = calc_interaction(p1_id, p2_id, 
                            system.dynamic_cfg, system)
                        
                        forces[p1_id] +=  f1
                        forces[p2_id] -=  f1
                    end
                end
            end
        end
    end
end

function calc_forces!(system::System, chunks::Chunks, device::Threaded)
    Threads.@threads :static for col in 1:chunks.num_cols
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
                    end
                end
            end
        end
    end

    # Soma todas as forças locais no array global
    get_forces(system) .= sum(system.forces)
end

"Compute total forces on particles."
function calc_forces!(system::System, chunks::Nothing, device::Sequencial)
    forces = get_forces(system)

    N = get_num_total_particles(system)
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

    for i in get_particles_ids(system)
        rel_pos = state.pos[i] - geometry_cfg.bottom_left

        out = (rel_pos .+ r) .> size_vec .|| (rel_pos .- r) .< 0
        if any(out)
            state.vel[i] = state.vel[i] .* (-2*out .+ 1)
        end
    end
end

function walls!(system::System, space_cfg::SpaceCfg{RigidWalls, G}) where G <: CircleCfg
    state = system.state
    circle = space_cfg.geometry_cfg
    max_r2 = (circle.radius - particle_radius(system.dynamic_cfg))^2
    for idx in get_particles_ids(system)
        pos = state.pos[idx]
        dr_2 = (pos - circle.center).^2
        r2 = sum(dr_2)
        if r2 <= max_r2
            continue
        end

        vel = state.vel[idx]
        new_vel = SVector(
            (vel.x*(dr_2.y - dr_2.x) - 2*vel.y*pos.x*pos.y) / r2,
            (-vel.y*(dr_2.y - dr_2.x) - 2*vel.x*pos.x*pos.y) / r2
        )
        state.vel[idx] = new_vel
    end
end

"Periodic walls."
function walls!(system::System, space_cfg::SpaceCfg{PeriodicWalls, G}) where G <: RectangleCfg
    state = system.state
    geometry_cfg = space_cfg.geometry_cfg
    center = geometry_cfg.bottom_left + geometry_cfg.size / 2
    half_size = geometry_cfg.size / 2

    for i in get_particles_ids(system)
        pos = state.pos[i]
        diff = pos - center
        
        out_bounds = abs.(diff) .> half_size
        if any(out_bounds)
            state.pos[i] = pos .- sign.(diff) .* (half_size * 2) .* out_bounds
        end
    end
end

"Slippery walls."
function walls!(system::System, space_cfg::SpaceCfg{SlipperyWalls, LinesCfg{T}}) where T
    state = system.state
    lines = space_cfg.geometry_cfg.lines
    p_radius = particle_radius(system.dynamic_cfg)
    
    for pid in get_particles_ids(system)
        pos_i = state.pos[pid] 
        
        for line in lines
            point = line.p1
            
            dr = pos_i - point
            delta_s = sum(dr .* line.normal)
            
            if abs(delta_s) > p_radius
                continue
            end
            
            delta_t = sum(dr .* line.tangent)
            
            is_corner = false
            if delta_t > 0
                if delta_t > line.length
                    if delta_t > (line.length + p_radius)
                        continue
                    end
                    is_corner = true
                    corner_p = line.p2
                end
            elseif delta_t > -p_radius    
                is_corner = true
                corner_p = line.p1
            else
                continue
            end
            
            if is_corner 
                dr = pos_i - corner_p
                norm = sqrt(sum(abs2, dr))
                if norm > p_radius
                    continue
                end
                alpha = p_radius / norm - 1
                state.pos[pid] += alpha * dr
            else
                sgn = sign(delta_s)
                alpha = sgn * (p_radius - sgn * delta_s) 
                state.pos[pid] += alpha * line.normal
            end
        end
    end
end

function walls!(system::System, space_cfg::SpaceCfg{SlipperyWalls, G}) where G <: CircleCfg
    pos = system.state.pos
    dynamic_cfg = system.dynamic_cfg
    state = system.state
    circle = space_cfg.geometry_cfg
    for pid in get_particles_ids(system)
        p = pos[pid]
        
        pr = get_particle_radius(dynamic_cfg, state, pid)
        max_r_2 = (circle.radius + pr)^2
        dr = calc_diff(p, circle.center, system.space_cfg)
        dr_2 = sum(dr.^2)
        if dr_2 > max_r_2
            continue
        end

        dr_norm = sqrt(dr_2)
        k = circle.radius + pr - dr_norm 

        pos[pid] = p + k * dr / dr_norm
    end
end

"Multiple Geometries"
function walls!(system::System, space_cfg::SpaceCfg{W, G}) where {W <: ManyWalls, G <: ManyGeometries}
    for (wall_cfg, geom_cfg) in zip(space_cfg.wall_type.list, space_cfg.geometry_cfg.list)
        walls!(system, SpaceCfg(wall_cfg, geom_cfg))
    end
end

"Quality of life."
@inline walls!(system::System) = walls!(system, system.space_cfg, system.type)
@inline walls!(system::System, space_cfg::SpaceCfg, ::Any) = walls!(system, space_cfg)

"Update state using Velocity-Verlet"
function update_verlet!(system::System, calc_forces_in!)
    state = system.state
    forces = get_forces(system)
    old_forces = copy(forces)

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

    for i in 1:get_num_total_particles(system)
        theta = state.pol_angle[i]
        pol = SVector(cos(theta), sin(theta) )

        vel = vo * pol + mu * forces[i]
        
        speed = sqrt(sum(abs, vel))
    
        cross_prod = (pol.x * vel.y - pol.y * vel.x) / speed 
        if abs(cross_prod) > 1
            cross_prod = sign(cross_prod)
        end

        d_theta = 1/relax_time * asin(cross_prod) * dt + sqrt(2 * dr * dt) * randn()
        
        state.pos[i] += vel * dt
        state.pol_angle[i] += d_theta
    end
end

function update_rtp!(system::System)
    state::SelfPropelledState = system.state
    forces = get_forces(system)

    dt = system.int_cfg.dt

    dynamic_cfg = system.dynamic_cfg
    vo, sigma, epsilon, tumble_rate = dynamic_cfg.vo, dynamic_cfg.sigma, dynamic_cfg.epsilon, dynamic_cfg.tumble_rate

    for i in 1:get_num_total_particles(system)
        theta = state.pol_angle[i]
        pol = SVector(cos(theta), sin(theta))

        vel = vo * pol + forces[i]
        # vel_x = vo * pol_x + forces[1, i]
        # vel_y = vo * pol_y + forces[2, i]
        
        # speed = sqrt(vel_x^2 + vel_y^2)
        speed = sqrt(sum(abs2, vel))
        
        # Update positions
        state.pos[i] += vel * dt
        # state.pos[1, i] += vel_x * dt
        # state.pos[2, i] += vel_y * dt

        # Update director
        u = rand()
        if u < tumble_rate * dt # accept tumble
            state.pol_angle[i] = 2*π*rand()
        end
    end
end

function update_time!(system)
    system.time_info.time += system.int_cfg.dt
    system.time_info.num_steps += 1
end


"Advance system one time step."
function newton_step!(system::System)
    clean_forces!(system)
    update_chunks!(system.chunks)
    calc_forces!(system)
    update_verlet!(system, calc_forces!)
    walls!(system, system.space_cfg)
    update_time!(system)
end

function szabo_step!(system::System)
    clean_forces!(system)
    update_chunks!(system.chunks)
    calc_forces!(system)
    update_szabo!(system)
    walls!(system, system.space_cfg)
    update_time!(system)
end

function rtp_step!(system::System)
    clean_forces!(system)
    update_chunks!(system.chunks)
    calc_forces!(system)
    update_rtp!(system)
    walls!(system, system.space_cfg)
    update_time!(system)
end

function get_step_function(::StandardSys, system)
    cfg_to_step = DefaultDict(
        newton_step!,
        SzaboCfg => szabo_step!,
        RunTumbleCfg => rtp_step!,
    )
    return cfg_to_step[typeof(system.dynamic_cfg)]
end

end