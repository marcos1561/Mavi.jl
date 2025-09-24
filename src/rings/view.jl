module RingsGraphs

using GLMakie

using Mavi.Rings
using Mavi.Rings.States
using Mavi.Rings.Configs: get_interaction_cfg, get_num_particles
import Mavi.Visualization.SystemGraphs: Graph, GraphCfg, MainGraph, MainGraphCfg, GraphComp, GraphCompCfg, get_graph_data, update_graph_data, get_graph, update_graph

function update_types_to_ring_id!(types, system)
    for ring_id in axes(system.rings_pos, 2)
        for p_id in 1:get_num_particles(system, ring_id)
            idx = to_scalar_idx(system.state, ring_id, p_id)
            types[idx] = ring_id
        end
    end
end

function get_graph_data(cfg::GraphCfg, system, state::RingsState)
    num_p = length(state.pos)

    types = Vector{Int}(undef, num_p)
    pos = similar(state.pos)

    idx = 1
    for ring_id in get_active_ids(state)
        for p_id in 1:get_num_particles(system.dynamic_cfg, state, ring_id)
            pos[idx] = state.rings_pos[p_id, ring_id] 
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
    num_total_particles = 0
    for ring_id in get_active_ids(state)
        num_p = get_num_particles(system.dynamic_cfg, state, ring_id)
        num_total_particles += num_p
        for p_id in 1:num_p
            pos[idx] = state.rings_pos[p_id, ring_id] 
            idx += 1
        end
    end

    cfg.pos_obs[] = @view pos[1:num_total_particles]
end

function get_graph_data(cfg::GraphCompCfg, system, state::RingsState)
    types = Vector{Int}(undef, length(state.pos))
    num_max_p = num_max_particles(state)
    for ring_id in axes(state.rings_pos, 2)
        for particle_id in 1:num_max_p
            types[to_scalar_idx(state, ring_id, particle_id)] = ring_id
        end
    end
    return types
end

@kwdef struct InvasionsGraphCfg{C, F} <: GraphCompCfg
    color::C = "black"
    update_data::F = update_graph_data
end

struct InvasionsGraph{O, C} <: GraphComp
    obs_list::O
    cfg::C
end

function get_graph(ax, pos_obs, system, cfg::InvasionsGraphCfg)
    invasions_pos_obs = Observable(eltype(system.state.pos)[])

    scatter!(ax, invasions_pos_obs, color=cfg.color, marker=:cross, markersize=12)
    graph = InvasionsGraph((invasions=invasions_pos_obs,), cfg)
    update_graph(graph, system)
    return graph
end

function update_graph_data(graph::InvasionsGraph, system)
    graph.obs_list[:invasions][] = [system.state.pos[inv.p_id] for inv in system.info.invasions.list]
    return
end


end # RingsGraphs