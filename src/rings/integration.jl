module Integration
    
import Mavi.Integration: calc_diff_and_dist, calc_interaction, calc_forces_chunks!
using Mavi.Rings.Configs
using Mavi.Rings.States
import Mavi.ChunksMod: Chunks, update_chunks!, update_particle_chunk!
using Mavi.Configs: SpaceCfg, RectangleCfg, PeriodicWalls

export step!

function to_scalar_idx(dynamic_cfg, ring_id, particle_id)
    return (ring_id - 1) * dynamic_cfg.num_max_particles + particle_id
end

function get_ring_id(idx, num_particles)
    div(idx-1, num_particles) + 1
end

function get_active_ids(state)
    state.active_ids[1:state.num_active]
end

function update_chunks!(chunks::Chunks{RingsState})
    for ring_id in get_active_ids(chunks.state)
        ring_type = chunks.state.type[ring_id]
        for p_id in 1:chunks.extra_info.num_particles[ring_type]
            update_particle_chunk!(chunks, to_scalar_idx(chunks.extra_info, ring_id, p_id))
        end
    end
end

function calc_interaction(i, j, state, dynamic_cfg::RingsCfg, space_cfg)
    r1 = get_ring_id(i, dynamic_cfg.num_max_particles)
    r2 = get_ring_id(j, dynamic_cfg.num_max_particles)

    t1, t2 = state.type[r1], state.type[r2]
    
    interaction_cfg = get_interaction_cfg(t1, t2, dynamic_cfg.interaction_finder)
    return calc_interaction_force(i, j, state, dynamic_cfg, interaction_cfg, space_cfg)
end

function calc_interaction_force(i, j, state, dynamic_cfg, interaction_cfg::HarmTruncCfg, space_cfg)
    dx, dy, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)
    
    if dist > interaction_cfg.dist_max
        return [0.0, 0.0]
    end

    num_max_particles = size(state.rings_pos)[2]
    ring_id1 = get_ring_id(i, num_max_particles)
    ring_id2 = get_ring_id(j, num_max_particles)
    
    if ring_id1 == ring_id2
        diff = abs(i - j)
        num_p = dynamic_cfg.num_particles[state.type[ring_id1]]
        
        # if diff % (state.num_particles[ring_id1] - 2) == 1 
        if diff == 1 || diff == num_p - 1
            return [0.0, 0.0]     
        end
    end

    dist_eq = interaction_cfg.dist_eq

    if dist < dist_eq
        # compute rep force
        fmod = -interaction_cfg.k_rep/dist_eq * (dist - dist_eq) # restoring force
        fx = fmod * dx / dist 
        fy = fmod * dy / dist
        
        return [fx, fy]
    end
    
    if ring_id1 == ring_id2
        return [0.0, 0.0]
    end
    
    # compute atr force
    adh_size = interaction_cfg.dist_max - dist_eq 
    fmod = -interaction_cfg.k_atr/adh_size*(dist - dist_eq) # restoring force
    fx = fmod * dx / dist # unit vector x_ij/dist
    fy = fmod * dy / dist
    
    return [fx, fy]
end

function springs_force(first_id, second_id, ring_id, state, rings_cfg, space_cfg)
    p1_id = to_scalar_idx(rings_cfg, ring_id, first_id)
    p2_id = to_scalar_idx(rings_cfg, ring_id, second_id)
    dx, dy, dist = calc_diff_and_dist(p1_id, p2_id, state.pos, space_cfg)
    
    ring_type = state.type[ring_id]
    k = rings_cfg.k_spring[ring_type]
    l = rings_cfg.l_spring[ring_type]

    fmod = -k * (dist - l)         

    spring_fx = dx/dist * fmod
    spring_fy = dy/dist * fmod

    return [spring_fx, spring_fy]
end

function cross_prod(v1, v2)
    return v1[1] * v2[2] - v1[2] * v2[1]
end

function calc_area(ring_pos)
    area = 0.0
    num_points = size(ring_pos)[2] 

    for i in 1:(num_points - 1)
        area += cross_prod(ring_pos[:, i], ring_pos[:, i+1])
    end
    area += cross_prod(ring_pos[:, end], ring_pos[:, 1])
    return area / 2.0
end

function update_continuos_pos!(system)
    continuos_pos = system.info.continuos_pos
    for ring_id in get_active_ids(system.state)
        continuos_pos[:, 1, ring_id] = system.state.rings_pos[:, 1, ring_id]
        ring_t = system.state.type[ring_id]
        for i in 2:system.dynamic_cfg.num_particles[ring_t]
            i_idx = to_scalar_idx(system.dynamic_cfg, ring_id, i)
            last_idx = to_scalar_idx(system.dynamic_cfg, ring_id, i-1)
            dx, dy, _ = calc_diff_and_dist(i_idx, last_idx, system.state.pos, system.space_cfg)
            
            continuos_pos[1, i, ring_id] = continuos_pos[1, i-1, ring_id] + dx
            continuos_pos[2, i, ring_id] = continuos_pos[2, i-1, ring_id] + dy
        end
    end
end

function get_continuos_pos(ring_id, system, wall_type)
    ring_type = system.state.type[ring_id]
    num_particles = system.dynamic_cfg.num_particles[ring_type]
    return @view system.state.ring_pos[:, 1:num_particles, ring_id]
end

function get_continuos_pos(ring_id, system, wall_type::PeriodicWalls)
    ring_type = system.state.type[ring_id]
    num_particles = system.dynamic_cfg.num_particles[ring_type]
    return @view system.info.continuos_pos[:, 1:num_particles, ring_id]
end

function area_forces!(system)
    rings_pos = system.state.rings_pos
    for ring_id in get_active_ids(system.state)
        area = calc_area(get_continuos_pos(ring_id, system, system.space_cfg.wall_type))
        system.info.areas[ring_id] = area

        ring_t = system.state.type[ring_id]
        k_area = system.dynamic_cfg.k_area[ring_t]
        p0 = system.dynamic_cfg.p0[ring_t]
        l0 = system.dynamic_cfg.l_spring[ring_t]
        num_particles = system.dynamic_cfg.num_particles[ring_t]

        area0 = (num_particles * l0 / p0)^2

        for i in 1:num_particles
            fmod = k_area * (area - area0)

            id1 = i - 1
            if id1 == 0
                id1 = num_particles
            end
            
            id2 = i + 1
            if id2 == num_particles + 1
                id2 = 1
            end
            
            id1_scalar = to_scalar_idx(system.dynamic_cfg, ring_id, id1)
            id2_scalar = to_scalar_idx(system.dynamic_cfg, ring_id, id2)
            dx, dy, _ = calc_diff_and_dist(id2_scalar, id1_scalar, system.state.pos, system.space_cfg)

            a_deriv_x = 1/2 * dy
            a_deriv_y = -1/2 * dx

            fx = -fmod * a_deriv_x
            fy = -fmod * a_deriv_y
            
            idx_scalar = to_scalar_idx(system.dynamic_cfg, ring_id, i)
            system.forces[1, idx_scalar] += fx
            system.forces[2, idx_scalar] += fy
        end
    end
end

function forces!(system)
    calc_forces_chunks!(system)

    state = system.state
    dynamic_cfg = system.dynamic_cfg
    space_cfg = system.space_cfg

    for ring_id in get_active_ids(state)
        ring_type = state.type[ring_id]
        num_particles = dynamic_cfg.num_particles[ring_type]
        for spring_id in 1:num_particles
            first_id = spring_id
            second_id = spring_id + 1
            if spring_id == num_particles
                second_id = 1
            end

            fx, fy = springs_force(first_id, second_id, ring_id, state, dynamic_cfg, space_cfg) 
            
            p1_id = to_scalar_idx(system.dynamic_cfg, ring_id, first_id)
            p2_id = to_scalar_idx(system.dynamic_cfg, ring_id, second_id)

            system.forces[1, p1_id] += fx
            system.forces[2, p1_id] += fy        
            system.forces[1, p2_id] -= fx
            system.forces[2, p2_id] -= fy        
        end
    end

    area_forces!(system)
end

"Periodic walls"
function walls!(system, space_cfg::SpaceCfg{PeriodicWalls, RectangleCfg}) 
    state = system.state
    geometry_cfg = space_cfg.geometry_cfg
    center = (geometry_cfg.bottom_left[1] + geometry_cfg.length/2, geometry_cfg.bottom_left[2] + geometry_cfg.height/2)

    for ring_id in get_active_ids(system.state)
        ring_type = system.state.type[ring_id]
        for i in 1:system.dynamic_cfg.num_particles[ring_type]
            pos_i = @view state.rings_pos[:, i, ring_id] 
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
end

function update!(system)
    state = system.state
    
    dt = system.int_cfg.dt

    state = system.state
    dynamic_cfg = system.dynamic_cfg
    
    for ring_id in get_active_ids(state)
        ring_type = state.type[ring_id]
        num_particles = dynamic_cfg.num_particles[ring_type]
        
        vo, relax_time = dynamic_cfg.vo[ring_type], dynamic_cfg.relax_time[ring_type]
        mu, dr = dynamic_cfg.mobility[ring_type], dynamic_cfg.rot_diff[ring_type] 
        
        theta = state.pol[ring_id]
        pol_x, pol_y = cos(theta), sin(theta) 
        vel_cm_x = 0.0
        vel_cm_y = 0.0

        for i in 1:num_particles

            p_id = to_scalar_idx(dynamic_cfg, ring_id, i)

            vel_x = vo * pol_x + mu * system.forces[1, p_id]
            vel_y = vo * pol_y + mu * system.forces[2, p_id]
            
            vel_cm_x += vel_x
            vel_cm_y += vel_y

            state.pos[1, p_id] += vel_x * dt
            state.pos[2, p_id] += vel_y * dt
        end
        vel_cm_x /= num_particles
        vel_cm_y /= num_particles

        speed = (vel_cm_x^2 + vel_cm_y^2)^.5

        cross_prod = (pol_x * vel_cm_y - pol_y * vel_cm_x) / speed 
        if abs(cross_prod) > 1
            cross_prod = sign(cross_prod)
        end
        d_theta = 1/relax_time * asin(cross_prod) * dt + sqrt(2 * dr * dt) * randn()
        state.pol[ring_id] += d_theta
    end
end

function step!(system, int_cfg::RingsIntCfg)
    update_chunks!(system.chunks)
    update_continuos_pos!(system)
    forces!(system)
    update!(system)
    walls!(system, system.space_cfg)
end

end