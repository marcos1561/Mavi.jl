using Mavi.Rings
using Mavi.Rings.States
using Mavi.Rings.Configs: get_interaction_cfg
import Mavi.Visualization.SystemGraphs: MainGraph, MainGraphCfg, Graph, GraphCfg, get_graph_data, update_graph_data

function get_graph_data(cfg::GraphCfg, system, state::RingsState)
    num_p = size(state.pos, 2)

    types = Vector{Int}(undef, num_p)
    pos = similar(state.pos)

    idx = 1
    for ring_id in get_active_ids(state)
        for p_id in 1:get_num_particles(system.dynamic_cfg, state, ring_id)
            pos[:, idx] .= state.rings_pos[:, p_id, ring_id] 
            types[idx] = ring_id
            idx += 1
        end
    end
    num_t = idx-1
    # pos_view = @view system.debug_info.graph_pos[:, 1:num_t]
    # types_view = @view system.debug_info.graph_type[1:num_t]

    return (pos=pos, types=types, num_p=num_t)
end

function update_graph_data(cfg::MainGraph, system, state::RingsState)
    pos = cfg.pos
    idx = 1
    for ring_id in get_active_ids(state)
        for p_id in 1:get_num_particles(system.dynamic_cfg, state, ring_id)
            pos[:, idx] .= state.rings_pos[:, p_id, ring_id] 
            idx += 1
        end
    end
end

function get_graph_data(system, state::RingsState)
    info = system.debug_info
    pos, radius, color = info.graph_pos, info.graph_radius, info.graph_color

    colors = [:red, :blue]

    idx = 1
    for ring_id in get_active_ids(state)
        ring_type = state.types[ring_id]
        for p_id in 1:system.dynamic_cfg.num_particles[ring_type]
            pos[1, idx] = state.rings_pos[1, p_id, ring_id] 
            pos[2, idx] = state.rings_pos[2, p_id, ring_id] 
            
            radius_i = get_interaction_cfg(ring_type, ring_type, system.dynamic_cfg.interaction_finder).dist_eq / 2.0
            radius[idx] = radius_i
            color[idx] = colors[ring_type]
            
            idx += 1
        end
    end

    count = idx-1
    return (pos=pos[:, 1:count], radius=radius[1:count], color=color[1:count])    
end