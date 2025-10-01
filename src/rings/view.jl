module RingsGraphs

export InvasionsGraphCfg, RingsNumsGraphCfg

using GLMakie

using Mavi.Rings
using Mavi.Rings.States
using Mavi.Rings.Configs
import Mavi.Visualization.SystemGraphs: Graph, GraphCfg, MainGraph, MainGraphCfg, GraphComp, GraphCompCfg, get_graph_data, update_graph_data, get_graph, update_graph

function update_types_to_ring_id!(types, system)
    for ring_id in axes(system.rings_pos, 2)
        for p_id in 1:get_num_particles(system, ring_id)
            idx = to_scalar_idx(system.state, ring_id, p_id)
            types[idx] = ring_id
        end
    end
end

function get_graph_data(cfg::GraphCfg, state::RingsState, system::Rings.System)
    num_p = length(state.pos)

    types = Vector{Int}(undef, num_p)
    pos = similar(state.pos)

    idx = 1
    for ring_id in get_rings_ids(state)
        for p_id in 1:ring_num_particles(state, ring_id)
            pos[idx] = state.rings_pos[p_id, ring_id] 
            types[idx] = ring_id
            idx += 1
        end
    end
    num_t = idx-1
    # pos_view = @view system.debug_info.graph_pos[:, 1:num_t]
    # types_view = @view system.debug_info.graph_type[1:num_t]

    return (pos=pos, types=types)
end

function update_graph_data(cfg::MainGraph, state::RingsState, system::Rings.System)
    pos = cfg.pos
    idx = 1
    num_total_particles = 0
    for ring_id in get_rings_ids(state)
        num_p = ring_num_particles(state, ring_id)
        num_total_particles += num_p
        for p_id in 1:num_p
            pos[idx] = state.rings_pos[p_id, ring_id] 
            idx += 1
        end
    end

    cfg.pos_obs[] = @view pos[1:num_total_particles]
end

function get_graph_data(cfg::GraphCompCfg, state::RingsState, system::Rings.System)
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

function update_graph_data(graph::InvasionsGraph, system::Rings.System)
    graph.obs_list[:invasions][] = [system.state.pos[inv.p_id] for inv in system.info.invasions.list]
    return
end

@kwdef struct RingsNumsGraphCfg <: GraphCompCfg 
    kwargs = Dict()
end

struct RingsNumsGraph{P} <: GraphComp 
    plot::P
    cfg::RingsNumsGraphCfg
end

function get_graph(ax, pos_obs, system, cfg::RingsNumsGraphCfg)
    kwargs = cfg.kwargs

    if !(:align in keys(kwargs))
        kwargs[:align] = (:center, :center)
    end

    plot = text!(ax, [zero(eltype(pos_obs[]))]; cfg.kwargs...)
    graph = RingsNumsGraph(plot, cfg)
    update_graph(graph, system)
    
    return graph
end

function update_graph(comp::RingsNumsGraph, system)
    rings_ids = get_rings_ids(system)
    
    points = Vector{Point2f}(undef, length(rings_ids))
    text = Vector{String}(undef, length(rings_ids))
    for (idx, rid) in enumerate(rings_ids)
        points[idx] = system.info.cms[rid]
        text[idx] = string(rid)
    end

    Makie.update!(comp.plot, points; text=text)
end


end # RingsGraphs