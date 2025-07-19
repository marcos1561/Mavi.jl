module Integration
    
import Mavi.Integration: calc_diff_and_dist, calc_interaction, calc_forces!, walls!
import Mavi.ChunksMod: Chunks, update_chunks!, update_particle_chunk!

using Mavi.Rings
using Mavi.Rings.NeighborsMod
using Mavi.Systems
using Mavi.Rings.Configs
using Mavi.Rings.States
using Mavi.Configs: SpaceCfg, RectangleCfg, LinesCfg, PeriodicWalls, SlipperyWalls

export step!
    
function update_chunks!(chunks::Chunks{T, StateT, InfoT}) where {T, StateT<:RingsState, InfoT}
    chunks.num_particles_in_chunk .= 0
    for ring_id in get_active_ids(chunks.state)
        for p_id in 1:get_num_particles(chunks.extra_info, chunks.state, ring_id)
            update_particle_chunk!(chunks, to_scalar_idx(chunks.state, ring_id, p_id))
        end
    end
end

function calc_num_neighbors(neighbors_count::Nothing, diff_dist, interaction_cfg, state) end

function calc_interaction(i, j, dynamic_cfg::RingsCfg, system::System)
    ri = get_ring_id(i, num_max_particles(system.state))
    rj = get_ring_id(j, num_max_particles(system.state))

    interaction_cfg = get_interaction_cfg(ri, rj, system.state, dynamic_cfg.interaction_finder)
    diff_dist = calc_diff_and_dist(i, j, system.state.pos, system.space_cfg)
    
    dist = diff_dist[end]
    max_dist = 2 * Configs.particle_radius(interaction_cfg)
    neigh_update!(system.info.p_neigh, i, j, dist, max_dist)
    neigh_update!(system.info.r_neigh, i, j, dist, max_dist)
    # calc_num_neighbors(system.info.neighbors_count, diff_dist, interaction_cfg, system.state) 

    return calc_interaction_force(i, j, ri, rj, diff_dist, interaction_cfg, system)
end

function calc_interaction_force(i, j, ring_id1, ring_id2, diff_dist, interaction_cfg::HarmTruncCfg, system::System)
    state = system.state
    dynamic_cfg = system.dynamic_cfg

    dx, dy, dist = diff_dist
    
    if dist > interaction_cfg.dist_max
        return 0.0, 0.0
    end

    # num_max_particles = num_max_particles(state)
    
    if ring_id1 == ring_id2
        diff = abs(i - j)
        num_p = get_num_particles(dynamic_cfg, state, ring_id1) 
        
        if diff == 1 || diff == num_p - 1
            return 0.0, 0.0     
        end
    end

    dist_eq = interaction_cfg.dist_eq

    if dist < dist_eq
        # compute rep force
        fmod = -interaction_cfg.k_rep/dist_eq * (dist - dist_eq) # restoring force
        fx = fmod * dx / dist 
        fy = fmod * dy / dist
        
        return fx, fy
    end
    
    if ring_id1 == ring_id2
        return 0.0, 0.0
    end
    
    # compute atr force
    adh_size = interaction_cfg.dist_max - dist_eq
    fmod = -interaction_cfg.k_atr/adh_size*(dist - dist_eq) # restoring force
    fx = fmod * dx / dist # unit vector x_ij/dist
    fy = fmod * dy / dist
    
    return fx, fy
end

function springs_force(first_id, second_id, ring_id, state, rings_cfg, space_cfg)
    p1_id = to_scalar_idx(state, ring_id, first_id)
    p2_id = to_scalar_idx(state, ring_id, second_id)
    dx, dy, dist = calc_diff_and_dist(p1_id, p2_id, state.pos, space_cfg)
    
    # k = rings_cfg.k_spring
    # l = rings_cfg.l_spring
    if !has_types_func(state)
        k = rings_cfg.k_spring
        l = rings_cfg.l_spring
    else
        ring_type = state.types[ring_id]
        k = rings_cfg.k_spring[ring_type]
        l = rings_cfg.l_spring[ring_type]
    end

    fmod = -k * (dist - l)         

    spring_fx = dx/dist * fmod
    spring_fy = dy/dist * fmod

    return spring_fx, spring_fy
end

@inline function cross_prod(v1, v2)
    return v1[1] * v2[2] - v1[2] * v2[1]
end

function calc_area(ring_pos)
    area = 0.0
    num_points = size(ring_pos, 2) 

    for i in 1:(num_points - 1)
        p1 = @view ring_pos[:, i]
        p2 = @view ring_pos[:, i+1]
        area += cross_prod(p1, p2)
    end
    p1 = @view ring_pos[:, end]
    p2 = @view ring_pos[:, 1]
    area += cross_prod(p1, p2)
    return area / 2.0
end

function update_continuos_pos!(system, wall_type) end

function update_continuos_pos!(system, wall_type::PeriodicWalls)
    continuos_pos = system.info.continuos_pos
    Threads.@threads for ring_id in get_active_ids(system.state)
        continuos_pos[:, 1, ring_id] = system.state.rings_pos[:, 1, ring_id]
        for i in 2:get_num_particles(system.dynamic_cfg, system.state, ring_id) 
            i_idx = to_scalar_idx(system.state, ring_id, i)
            last_idx = to_scalar_idx(system.state, ring_id, i-1)
            dx, dy, _ = calc_diff_and_dist(i_idx, last_idx, system.state.pos, system.space_cfg)
            
            continuos_pos[1, i, ring_id] = continuos_pos[1, i-1, ring_id] + dx
            continuos_pos[2, i, ring_id] = continuos_pos[2, i-1, ring_id] + dy
        end
    end
end

function get_continuos_pos(ring_id, system, wall_type)
    num_particles = get_num_particles(system.dynamic_cfg, system.state, ring_id) 
    return @view system.state.rings_pos[:, 1:num_particles, ring_id]
end

function get_continuos_pos(ring_id, system, wall_type::PeriodicWalls)
    num_particles = get_num_particles(system.dynamic_cfg, system.state, ring_id) 
    return @view system.info.continuos_pos[:, 1:num_particles, ring_id]
end

function area_forces!(system)
    has_types = has_types_func(system.state)
    if !has_types
        k_area = system.dynamic_cfg.k_area
        p0 = system.dynamic_cfg.p0
        l0 = system.dynamic_cfg.l_spring
        num_particles = get_num_particles(system.dynamic_cfg)
    end
    
    # Threads.@threads for ring_id in get_active_ids(system.state)
    #     area = calc_area(get_continuos_pos(ring_id, system, system.space_cfg.wall_type))
    #     system.info.areas[ring_id] = area
    # end

    forces = get_forces(system)

    for ring_id in get_active_ids(system.state)
        area = calc_area(get_continuos_pos(ring_id, system, system.space_cfg.wall_type))
        system.info.areas[ring_id] = area
        # area = system.info.areas[ring_id]

        if has_types
            ring_type = system.state.types[ring_id]
            k_area = system.dynamic_cfg.k_area[ring_type]
            p0 = system.dynamic_cfg.p0[ring_type]
            l0 = system.dynamic_cfg.l_spring[ring_type]
            num_particles = system.dynamic_cfg.num_particles[ring_type]
        end

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
            
            id1_scalar = to_scalar_idx(system.state, ring_id, id1)
            id2_scalar = to_scalar_idx(system.state, ring_id, id2)
            @inbounds dx, dy, _ = calc_diff_and_dist(id2_scalar, id1_scalar, system.state.pos, system.space_cfg)

            a_deriv_x = 1/2 * dy
            a_deriv_y = -1/2 * dx

            fx = -fmod * a_deriv_x
            fy = -fmod * a_deriv_y
            
            idx_scalar = to_scalar_idx(system.state, ring_id, i)
            forces[1, idx_scalar] += fx
            forces[2, idx_scalar] += fy
        end
    end
end

function forces!(system)
    calc_forces!(system) 

    state = system.state
    dynamic_cfg = system.dynamic_cfg
    space_cfg = system.space_cfg
    forces = get_forces(system)

    for ring_id in get_active_ids(state)
        num_particles = get_num_particles(dynamic_cfg, state, ring_id)
        for spring_id in 1:num_particles
            first_id = spring_id
            second_id = spring_id + 1
            if spring_id == num_particles
                second_id = 1
            end

            fx, fy = springs_force(first_id, second_id, ring_id, state, dynamic_cfg, space_cfg) 
            
            p1_id = to_scalar_idx(system.state, ring_id, first_id)
            p2_id = to_scalar_idx(system.state, ring_id, second_id)

            forces[1, p1_id] += fx
            forces[2, p1_id] += fy        
            forces[1, p2_id] -= fx
            forces[2, p2_id] -= fy        
        end
    end

    area_forces!(system)
end

"Periodic walls"
function walls!(system::System, space_cfg::SpaceCfg{PeriodicWalls, RectangleCfg{T}}, dynamic_cfg::RingsCfg) where T 
    state = system.state
    geometry_cfg = space_cfg.geometry_cfg
    center = (geometry_cfg.bottom_left[1] + geometry_cfg.length/2, geometry_cfg.bottom_left[2] + geometry_cfg.height/2)

    for ring_id in get_active_ids(system.state)
        for i in 1:get_num_particles(system.dynamic_cfg, state, ring_id)
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

function walls!(system::System, space_cfg::SpaceCfg{SlipperyWalls, LinesCfg{T}}, dynamic_cfg::RingsCfg) where T
    state = system.state
    lines = space_cfg.geometry_cfg.lines
    p_radius = Configs.particle_radius(system.dynamic_cfg)
    for ring_id in get_active_ids(system.state)
        for i in 1:get_num_particles(system.dynamic_cfg, state, ring_id) 
            pos_i = @view state.rings_pos[:, i, ring_id] 
            
            for line in lines
                point_x = line.p1.x #+ line.normal.x * p_radius
                point_y = line.p1.y #+ line.normal.y * p_radius
                
                x, y = pos_i[1] - point_x, pos_i[2] - point_y
                delta_s = x * line.normal.x + y * line.normal.y   
                
                if abs(delta_s) > p_radius
                    continue
                end
                
                delta_t = x * line.tangent.x + y * line.tangent.y   
                
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
                    dx = pos_i[1] - corner_p.x
                    dy = pos_i[2] - corner_p.y
                    norm = sqrt(dx^2 + dy^2)
                    if norm > p_radius
                        continue
                    end
                    alpha = p_radius / sqrt(dx^2 + dy^2) - 1
                    pos_i[1] += alpha * dx
                    pos_i[2] += alpha * dy
                else
                    sgn = sign(delta_s)
                    alpha = sgn * (p_radius - sgn * delta_s) 
                    pos_i[1] += alpha * line.normal.x
                    pos_i[2] += alpha * line.normal.y
                end
            end
        end
    end
end

function update!(system)
    state = system.state
    forces = get_forces(system)

    dt = system.int_cfg.dt

    state = system.state
    dynamic_cfg = system.dynamic_cfg
    
    has_types = has_types_func(state) 

    if !has_types
        vo, relax_time = dynamic_cfg.vo, dynamic_cfg.relax_time
        mu, dr = dynamic_cfg.mobility, dynamic_cfg.rot_diff 
    end

    for ring_id in get_active_ids(state)
        num_particles = get_num_particles(dynamic_cfg, state, ring_id)
        
        if has_types
            ring_type = state.types[ring_id]
            vo, relax_time = dynamic_cfg.vo[ring_type], dynamic_cfg.relax_time[ring_type]
            mu, dr = dynamic_cfg.mobility[ring_type], dynamic_cfg.rot_diff[ring_type] 
        end

        theta = state.pol[ring_id]
        pol_x, pol_y = cos(theta), sin(theta) 
        vel_cm_x = 0.0
        vel_cm_y = 0.0

        for i in 1:num_particles
            p_id = to_scalar_idx(state, ring_id, i)

            vel_x = vo * pol_x + mu * forces[1, p_id]
            vel_y = vo * pol_y + mu * forces[2, p_id]
            
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

function cleaning!(system)
    clean_forces!(system)
    # neigh_clean!(system.info.p_neigh)
    # neigh_clean!(system.info.r_neigh)
end

function step!(system)
    update_chunks!(system.chunks)
    update_continuos_pos!(system, system.space_cfg.wall_type)
    
    cleaning!(system)
    forces!(system)
    # neigh_sum_buffers(system.info.p_neigh)
    # neigh_sum_buffers(system.info.r_neigh)
    
    update!(system)
    walls!(system)

    system.time_info.num_steps += 1
    system.time_info.time += system.int_cfg.dt
end

end