module SystemGraphs

export MainGraph, GraphCfg, GraphComp, GraphCompCfg, GraphCompDebug
export MainGraphCfg, CircleGraphCfg, ScatterGraphCfg
export drawn_borders, colors_from_cmap
export RingsGraphs

using GLMakie, ColorSchemes, DataStructures, Random
using Mavi.Systems
using Mavi.States
using Mavi.Configs

"Return circle vertices centered in the origin with given radius."
function get_circle(radius, resolution=20)
    theta = LinRange(0, 2π, resolution+1)

    circle = Vector{Point2f}()
    for id in 1:length(theta)-1
        t1, t2 = theta[id], theta[id+1]
        push!(circle, Point2f(radius * cos(t1), radius * sin(t1))) 
        push!(circle, Point2f(radius * cos(t2), radius * sin(t2))) 
    end

    return circle
end

function drawn_borders(ax, geometry_cfg::RectangleCfg; adjust_lims=true, color=:black)
    l, h = geometry_cfg.length, geometry_cfg.height
    x, y = geometry_cfg.bottom_left
    # lines!(ax, [x, y], [x+l, y], color=:black)
    # lines!(ax, [x+l, y], [x+l, y+h], color=:black)
    # lines!(ax, [x+l, y+h], [x, y+h], color=:black)
    # lines!(ax, [y+h, x], [x, y], color=:black)
    
    lines!(ax, [x, x+l], [y, y], color=color)
    lines!(ax, [x+l, x+l], [y+h, y], color=color)
    lines!(ax, [x+l, x], [y+h, y+h], color=color)
    lines!(ax, [x, x], [y+h, y], color=color)

    if adjust_lims
        dx = l*0.1
        dy = h*0.1
        xlims!(ax, x - dx, x + l + dx)
        ylims!(ax, y - dy, y + h + dy)
    end
end

function drawn_borders(ax, geometry_cfg::LinesCfg; color=:black)
    for line in geometry_cfg.lines
        p1, p2 = line.p1, line.p2
        lines!(ax, [p1.x, p2.x], [p1.y, p2.y], color=color)
    end
end

function drawn_borders(ax, geometry_cfg::CircleCfg)
    arc!(ax, geometry_cfg.center, geometry_cfg.radius, -π, π, color=:black)
end

function drawn_borders(ax, geometry_cfg::ManyGeometries)
    for g in geometry_cfg.list
        drawn_borders(ax, g)
    end
end

function update_circles!(circle_points, circle_base, pos)
    num_edges = length(circle_base)
    num_p = size(pos)[2]
    for i in 1:num_p
        begin_id = num_edges * (i-1) + 1
        end_id = begin_id + num_edges - 1 
        center = [pos[1, i], pos[2, i]]
        circle_points[begin_id:end_id] .= circle_base .+ [center]
    end
end

function get_num_types(types, cmap)
    if typeof(cmap) <: AbstractArray
        return max(length(types), length(cmap))
    else
        return length(types)
    end
end

function colors_from_cmap(cmap, num)
    if cmap isa Symbol
        cmap = colorschemes[cmap]
    end
    return Vector{Any}([cmap[rand(Float64)] for _ in 1:num])
end

function get_color_map(cmap::AbstractVector, num_types; rng=nothing) 
    if isnothing(rng)
        rng = Random.GLOBAL_RNG
    end

    cmap_out = Vector{RGBf}(undef, num_types)
    for (idx, c) in enumerate(cmap)
        if c isa RGBf
            c_rgb = c
        elseif c == :random
            c_rgb = rand(rng, RGBf)
        else
            c_rgb = RGBf(GLMakie.to_color(c))
        end
        cmap_out[idx] = c_rgb
    end

    for idx in (length(cmap)+1):num_types
       cmap_out[idx] = RGBf(GLMakie.to_color(:purple))
    end

    return cmap_out
end

function get_color_map(cmap::Symbol, num_types; rng=nothing)
    if isnothing(rng)
        rng = Random.GLOBAL_RNG
    end

    try
        cmap = colorschemes[cmap]
    catch end

    if cmap == :random
        return [rand(rng, RGBf) for _ in 1:num_types]
    elseif cmap == :ids
        return [rand(rng, RGBf) for _ in 1:num_types]
    elseif cmap isa ColorScheme
        return [cmap[rand(Float32)] for _ in 1:num_types]            
    else
        color = RGBf(GLMakie.to_color(cmap))
        return fill(color, num_types)
    end
end

function get_colors!(colors, types, cmap, particles_ids)
    count = 1
    for idx in particles_ids
        colors[count] = cmap[types[idx]]
        count += 1
    end

    cs = @view colors[1:count-1]
    return cs
end

function get_radius!(radius, dynamic_cfg, state, particles_ids)
    count = 1
    for idx in particles_ids
        radius[count] = get_particle_radius(dynamic_cfg, state, idx)
        count += 1
    end
    rs = @view radius[1:count-1]
    return rs
end


abstract type GraphCfg end
abstract type Graph end

function get_graph_data(graph_cfg::GraphCfg, state::State, system)
    num_p = length(system.state.pos)
    return (pos=state.pos, types=Vector(1:num_p))
end
@inline get_graph_data(graph::GraphCfg, system) = get_graph_data(graph, system.state, system)

update_graph_data(graph::Graph, state::State, system::System) = DefaultDict{Symbol, Bool}(false)
@inline update_graph_data(graph::Graph, system) = update_graph_data(graph, system.state, system)


abstract type GraphCompCfg <: GraphCfg end
abstract type GraphComp <: Graph end
abstract type GraphCompDebug <: GraphComp end

get_graph_data(graph_cfg::GraphCompCfg, state::State, system) = Vector(1:length(state.pos))

get_comp_update_data(c::GraphComp) = c.cfg.update_data
get_comp_update_data(c::GraphCompDebug) = (_, _) -> DefaultDict{Symbol, Bool}(false)

get_comp_obs_list(c::GraphComp) = c.obs_list
get_comp_obs_list(c::GraphCompDebug) = ()

function notify_comp_observables(comp::GraphComp, update_list)
    update_all = isnothing(update_list)
    for (name, obs) in pairs(get_comp_obs_list(comp))
        if update_all || update_list[name]
            notify(obs)
        end
    end
end

function update_graph(comp::GraphComp, system)
    update_list = get_comp_update_data(comp)(comp, system)
    notify_comp_observables(comp, update_list)
end


struct ScatterGraphCfg{C, KW, F} <: GraphCompCfg
    colors_map::C
    kwargs::KW
    update_data::F
    rng::AbstractRNG
end
function ScatterGraphCfg(;colors_map=:random, update_data=nothing, kwargs=nothing, rng=nothing)
    if colors_map isa Vector{Symbol}
        colors_map = RGBf.(GLMakie.to_color.(colors_map))
    end

    if update_data === nothing
        update_data = update_graph_data
    end

    if kwargs === nothing
        kwargs = ()
    end

    if isnothing(rng)
        rng = Random.GLOBAL_RNG
    end

    ScatterGraphCfg(colors_map, kwargs, update_data, rng)
end

struct ScatterGraph{O, C, P, PosObs} <: GraphComp
    types::Vector{Int}
    colors::Vector{C}
    cmap::Vector{C}
    plot::P
    obs_list::O
    pos_obs::PosObs
    cfg::ScatterGraphCfg
end

function get_graph(ax, pos_obs, system, cfg::ScatterGraphCfg)
    num_total_particles = length(system.state.pos)
    
    types = get_graph_data(cfg, system.state, system)
    cmap = get_color_map(cfg.colors_map, get_num_types(types, cfg.colors_map), rng=cfg.rng)
    colors = Vector{eltype(cmap)}(undef, num_total_particles)

    scatter_plot = scatter!(ax, [zero(eltype(pos_obs[]))]; cfg.kwargs...)

    comp = ScatterGraph(types, colors, cmap, scatter_plot, (), pos_obs, cfg)
    update_graph(comp, system)
    return comp
end

function update_graph(comp::ScatterGraph, system)
    update_list = get_comp_update_data(comp)(comp, system)
    notify_comp_observables(comp, update_list)

    particles_ids = get_particles_ids(system)
    colors = get_colors!(comp.colors, comp.types, comp.cmap, particles_ids)
    pos = comp.pos_obs[]
    points = [Point2f(p) for p in pos]

    Makie.update!(comp.plot, points; color=colors)
end

@kwdef struct StrokeCfg{C}
    color::C = "black"
    width::Float64 = 0.5
end

struct CircleGraphCfg{C, S<:Union{StrokeCfg, Nothing}, F} <: GraphCompCfg
    circle_radius::Float64
    stroke_cfg::S
    circle_rel::Int
    colors_map::C
    update_data::F
    rng::AbstractRNG
end
function CircleGraphCfg(;circle_radius=-1.0, stroke_cfg=StrokeCfg(), circle_rel=20, colors_map=:random, 
    update_data=nothing, rng=nothing)
    
    # if colors_map === nothing
    #     colors_map = RGBf(GLMakie.to_color(:black))
    # end

    if colors_map isa Vector{Symbol}
        colors_map = RGBf.(GLMakie.to_color.(colors_map))
    end

    if update_data === nothing
        update_data = update_graph_data
    end

    if isnothing(rng)
        rng = Random.GLOBAL_RNG
    end

    CircleGraphCfg(circle_radius, stroke_cfg, circle_rel, colors_map, update_data, rng)
end

struct CircleGraph{C, P, O, PosObs} <: GraphComp
    types::Vector{Int}
    colors::Vector{C}
    radius::Vector{Float64}
    cmap::Vector{C}
    plot::P
    obs_list::O
    pos_obs::PosObs
    cfg::CircleGraphCfg
end

function get_graph(ax, pos_obs, system, cfg::CircleGraphCfg)
    num_total_particles = length(system.state.pos)

    # TODO: Make cfg.circle_radius work with dispatch system
    if cfg.circle_radius == -1
        radius = Vector{eltype(system.state.pos[1])}(undef, num_total_particles)
    else
        if cfg.circle_radius <: Number
            radius = fill(cfg.circle_radius, num_total_particles) 
        elseif length(cfg.circle_radius) != num_total_particles
            throw(ArgumentError(
                "Length of `cfg.circle_radius` ($(length(cfg.circle_radius))) should be " *
                "the same as number of particles ($(num_total_particles))"
            ))
        end
    end

    types = get_graph_data(cfg, system.state, system)
    cmap = get_color_map(cfg.colors_map, get_num_types(types, cfg.colors_map), rng=cfg.rng)
    colors = Vector{eltype(cmap)}(undef, num_total_particles)

    if isnothing(cfg.stroke_cfg)
        circles_plot = poly!(ax, [Circle(Point2f(0, 0), 1)])
    else
        circles_plot = poly!(ax,
            [Circle(Point2f(0, 0), 1)],
            strokecolor=cfg.stroke_cfg.color,
            strokewidth=cfg.stroke_cfg.width,
        )
    end
    
    comp = CircleGraph(types, colors, radius, cmap, circles_plot, (), pos_obs, cfg)
    update_graph(comp, system)
    return comp
end

function update_graph(comp::CircleGraph, system)
    update_list = get_comp_update_data(comp)(comp, system)
    notify_comp_observables(comp, update_list)

    particles_ids = get_particles_ids(system)
    colors = get_colors!(comp.colors, comp.types, comp.cmap, particles_ids)
    radius = get_radius!(comp.radius, system.dynamic_cfg, system.state, particles_ids)

    pos = comp.pos_obs[]
    
    circles = [Circle(Point2f(p), r) for (p, r) in zip(pos, radius)]
    Makie.update!(comp.plot, circles; color=colors)
end

"""
Drawn every particle with a circumference and a central point.

- circle_radius:  
    Radius used to draw circles centered on particles positions.  
    If not given, `particle_radius(dynamic_cfg)` will be used.
- circle_rel:  
    How many vertices used to draw circles.
"""
struct MainGraphCfg{T<:Tuple} <: GraphCfg
    comps_cfgs::T
end
function MainGraphCfg(comps::AbstractVector)
    MainGraphCfg(tuple(comps...))
end
function MainGraphCfg(comp::GraphCompCfg)
    MainGraphCfg((comp,))
end
function MainGraphCfg()
    MainGraphCfg(CircleGraphCfg())
end

struct MainGraph{A, P, O, C<:Tuple} <: Graph
    ax::A
    pos::P
    pos_obs::O
    components::C
    cfg::MainGraphCfg
end

"Build the Graph respective to `cfg` using the given `grid_layout` from Makie."
function get_graph(ax, system, cfg::MainGraphCfg)
    data = get_graph_data(cfg, system.state, system)
    pos_obs = Observable(data.pos)

    drawn_borders(ax, system.space_cfg.geometry_cfg)
    
    graphs = Vector{GraphComp}()
    for g in cfg.comps_cfgs
        push!(graphs, get_graph(ax, pos_obs, system, g))
    end

    graph = MainGraph(ax, similar(data.pos), pos_obs, tuple(graphs...), cfg)
    update_graph(graph, system)
    graph
end

function update_graph_data(graph::MainGraph, part_ids, state::State) end
function update_graph_data(graph::MainGraph, part_ids::States.AbstractParticleIds, state::State)
    # pos = cfg.pos
    # idx = 1
    # num_total_particles = 0
    # for ring_id in get_rings_ids(state)
    #     num_p = get_num_particles(system.dynamic_cfg, state, ring_id)
    #     num_total_particles += num_p
    #     for p_id in 1:num_p
    #         pos[idx] = state.rings_pos[p_id, ring_id] 
    #         idx += 1
    #     end
    # end

    pos = graph.pos
    for (idx, p_id) in enumerate(get_particles_ids(state))
        pos[idx] = state.pos[p_id]
    end

    graph.pos_obs[] = @view pos[1:get_num_total_particles(state)]
end

function update_graph_data(graph::MainGraph, state::State, system::System)
    update_graph_data(graph, state.part_ids, state)
end

function update_graph(graph::MainGraph, system)
    update_graph_data(graph, system)
    for c in graph.components
        update_graph(c, system)
    end
    notify(graph.pos_obs)
end

include("../rings/view.jl")

end