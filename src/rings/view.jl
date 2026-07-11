module RingsGraphs

export InvasionsGraphCfg, RingsNumsGraphCfg, RingForceGraphCfg, NeighborsGraphCfg

using GLMakie

using Mavi.Rings
using Mavi.Rings.States
using Mavi.Rings.Configs
using Mavi.Rings.RingsDebug
using Mavi.Rings.NeighborsMod
import Mavi.Visualization.SystemGraphs: 
    Graph, GraphCfg, MainGraph, MainGraphCfg, 
    GraphComp, GraphCompCfg, 
    get_graph_data, update_graph_data, 
    get_graph, update_graph, 
    get_default_num_types

function update_types_to_ring_id!(types, system)
    for ring_id in axes(system.rings_pos, 2)
        for p_id in 1:get_num_particles(system, ring_id)
            idx = to_scalar_idx(system.state, ring_id, p_id)
            types[idx] = ring_id
        end
    end
end

function get_graph_data(cfg::GraphCfg, state::RingsState)
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

function update_graph_data(cfg::MainGraph, state::RingsState)
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

function get_graph_data(cfg::GraphCompCfg, state::RingsState)
    types = Vector{Int}(undef, length(state.pos))
    num_max_p = num_max_particles(state)
    for ring_id in axes(state.rings_pos, 2)
        for particle_id in 1:num_max_p
            types[to_scalar_idx(state, ring_id, particle_id)] = ring_id
        end
    end
    return types
end

function get_default_num_types(cfg::GraphCompCfg, state::RingsState)
    return size(state.rings_pos, 2)
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


# = 
# Single Ring Forces
# =
struct ForceCfg{C}
    color::C
    lengthscale::Float64
    kwargs::Dict
end
function ForceCfg(name; color=nothing, lengthscale=nothing, kwargs=())
    kwargs = Dict(kwargs)

    if lengthscale === nothing
        lengthscale = 1.0
    end

    name_to_color = (
        area="green",
        spring="red",
        pairwise="blue",
    )
    if color === nothing
        color = name_to_color[name]
    end

    kwargs = merge(Dict(kwargs), Dict(:color=>color, :lengthscale=>lengthscale))
    ForceCfg(color, lengthscale, kwargs)
end

struct RingForceGraphCfg <: GraphCfg
    ring_forces::DebugRingForces
    area_cfg::ForceCfg
    spring_cfg::ForceCfg
    pairwise_cfg::ForceCfg
    show_list::Tuple
    pairwise_mode::Symbol
end
function RingForceGraphCfg(;
    ring_forces,
    area_cfg=nothing,
    spring_cfg=nothing,
    pairwise_cfg=nothing,
    show_list=(:area, :spring, :pairwise),
    lengthscale=nothing,
    pairwise_mode=:all,
)
    if area_cfg === nothing
        area_cfg = ForceCfg(:area, lengthscale=lengthscale)
    end
    if spring_cfg === nothing
        spring_cfg = ForceCfg(:spring, lengthscale=lengthscale)
    end
    if pairwise_cfg === nothing
        pairwise_cfg = ForceCfg(:pairwise, lengthscale=lengthscale)
    end
    RingForceGraphCfg(ring_forces, area_cfg, spring_cfg, pairwise_cfg, show_list, pairwise_mode)
end

struct RingForceGraph{P, C} <: GraphComp
    arrows::P
    cfg::C
end

function get_graph(ax, system, cfg::RingForceGraphCfg)
    plot_pairwise = arrows2d!(ax, [0.0], [0.0], [0.0], [0.0]; cfg.pairwise_cfg.kwargs...)
    plot_spring = arrows2d!(ax, [0.0], [0.0], [0.0], [0.0]; cfg.spring_cfg.kwargs...)
    plot_area = arrows2d!(ax, [0.0], [0.0], [0.0], [0.0]; cfg.area_cfg.kwargs...)
    plots = (area=plot_area, spring=plot_spring, pairwise=plot_pairwise)
    graph = RingForceGraph(plots, cfg)
    update_graph(graph, system)
    return graph
end

function update_graph(graph::RingForceGraph, system)
    pos = system.state.pos
    ring_forces = graph.cfg.ring_forces
    ring_id = ring_forces.ring_id
    n = States.ring_num_particles(system, ring_id)

    pairwise_forces = zeros(eltype(ring_forces.area_forces), n)
    if graph.cfg.pairwise_mode === :all
        for forces_from_another_ring in ring_forces.interaction_forces
            pairwise_forces .+= forces_from_another_ring[1:n]
        end
    end

    forces = (spring=ring_forces.spring_forces, area=ring_forces.area_forces, pairwise=pairwise_forces)

    for name in (:spring, :area, :pairwise)
        x = Vector{Float64}(undef, n)
        y = Vector{Float64}(undef, n)
        u = Vector{Float64}(undef, n)
        v = Vector{Float64}(undef, n)

        forces_i = forces[name]
        for p_i in 1:n
            idx = States.to_scalar_idx(system.state, ring_id, p_i)
            x[p_i] = pos[idx][1]
            y[p_i] = pos[idx][2]
            u[p_i] = forces_i[p_i][1]
            v[p_i] = forces_i[p_i][2]
        end
        Makie.update!(graph.arrows[name], x, y, u, v)
    end
end

# ==
# Neighbors
# ==
@kwdef struct NeighborsGraphCfg{N} <: GraphCfg 
    neighbors::N
    offset=nothing
end

struct NeighborsGraph{N, P} <: Graph 
    cfg::NeighborsGraphCfg{N}
    plot::P
end

function get_graph(ax, system, cfg::NeighborsGraphCfg)
    plot = text!(ax, [zero(eltype(system.state.pos))])
    graph = NeighborsGraph(cfg, plot)
    update_graph(graph, system)
    graph
end

function update_graph(graph::NeighborsGraph, system) 
    neigh = get_neigh(graph.cfg.neighbors)

    pos = system.state.pos
    points = Point2f[]
    text = String[]
    count = get_neigh_count(neigh)
    
    offset = graph.cfg.offset
    if offset === nothing
        offset = zero(eltype(pos))
    end
        
    for (pid, c) in enumerate(count)
        push!(points, pos[pid] + offset)
        push!(text, "$c")
    end

    Makie.update!(graph.plot, points; text=text)
end

end # RingsGraphs