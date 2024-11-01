module SystemGraphs

export GraphCfg, DefaultGraphCfg
export CircleGraph, CircleGraphCfg

using GLMakie
using Mavi: State, System
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

function drawn_borders(ax, space_cfg::RectangleCfg)
    l, h = space_cfg.length, space_cfg.height
    lines!(ax, [0, l], [0, 0], color=:black)
    lines!(ax, [0, l], [h, h], color=:black)
    lines!(ax, [0, 0], [0, h], color=:black)
    lines!(ax, [l, l], [0, h], color=:black)

    xlims!(ax, -l*0.1, l*1.1)
    ylims!(ax, -h*0.1, h*1.1)
end

function drawn_borders(ax, space_cfg::CircleCfg)
    arc!(ax, Point2f(0), space_cfg.radius, -π, π, color=:black)
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


abstract type GraphCfg end

"""
Drawn every particle with a circumference and a central point.

- circle_radius:  
    Radius used to draw circles centered on particles positions.  
    If not given, `particle_radius(dynamic_cfg)` will be used.
- circle_rel:  
    How many vertices used to draw circles.
"""
@kwdef struct DefaultGraphCfg <: GraphCfg 
    circle_radius = -1.0
    circle_rel = 20
    colors_map = [:blue]
end

struct DefaultGraph{T1, T2, T3}
    ax::T1
    circle_base::Vector{Point2f}
    circle_points::Vector{Point2f}
    circle_points_obs::Observable{Vector{Point2f}}
    pos_obs::T2
    colors::Vector{Symbol}
    colors_obs::T3
    cfg::DefaultGraphCfg
end

function get_graph_data(graph_cfg::DefaultGraphCfg, system, state)
    return (pos=state.pos, types=fill(1, size(state.pos)[2]))
end

function get_colors(graph_cfg::DefaultGraphCfg, colors, types)
    count = 1
    for t in types
        colors[count] = graph_cfg.colors_map[t]
        count += 1
    end

    cv = @view colors[1:length(types)]
    return cv
end

"Build the Graph respective to `cfg` using the given `grid_layout` from Makie."
function get_graph(grid_layout, system, cfg::DefaultGraphCfg)
    ax = Axis(grid_layout[1, 1])

    data = get_graph_data(cfg, system, system.state)
    pos = data.pos

    x_obs = Observable(pos[1, :])
    y_obs = Observable(pos[2, :])
    pos_obs = [x_obs, y_obs]

    colors = Vector{Symbol}(undef, size(pos)[2])
    colors_obs = Observable(get_colors(cfg, colors, data.types))

    scatter!(ax, x_obs, y_obs, color=colors_obs)
    # ax.aspect = 1

    drawn_borders(ax, system.space_cfg.geometry_cfg)

    if cfg.circle_radius == -1
        radius = particle_radius(system.dynamic_cfg) 
    else
        radius = cfg.circle_radius 
    end

    # Circles
    num_p = size(system.state.pos)[2]
    circle_base = get_circle(radius, cfg.circle_rel)
    circle_points = Vector{Point2f}(undef, num_p * length(circle_base))

    update_circles!(circle_points, circle_base, pos)
    circle_points_obs = Observable(circle_points)
    linesegments!(ax, circle_points_obs, color=:black)

    DefaultGraph(ax, circle_base, circle_points, circle_points_obs, pos_obs,
        colors, colors_obs, cfg,
    )
end

function update_graph(graph::DefaultGraph, system)
    data = get_graph_data(graph.cfg, system, system.state)
    num_p = size(data.pos)[2]
    
    update_circles!(graph.circle_points, graph.circle_base, data.pos)
    graph.circle_points_obs[] = graph.circle_points[1:num_p * length(graph.circle_base)]
    
    graph.pos_obs[1].val = data.pos[1, :]
    graph.pos_obs[2].val = data.pos[2, :]
    graph.colors_obs.val = get_colors(graph.cfg, graph.colors, data.types)
    
    notify(graph.pos_obs[1])
    # notify(graph.pos_obs[2])
    
    # notify(graph.circle_points_obs)
    # notify(graph.pos_obs[1])
    # notify(graph.pos_obs[2])
end

"""
Drawn every particle as circle with given radius and color.
"""
struct CircleGraphCfg <: GraphCfg end

mutable struct CircleGraph{T1, T2}
    ax::T1
    circles::T2
end

function get_graph_data(system, state) end

function update_circles!(graph, pos, radius, color)
    for c in graph.circles
        delete!(graph.ax, c)
    end

    graph.circles = [mesh!(graph.ax, Circle(p, r), color=c, shading=NoShading) for (p, r, c) in zip(pos, radius, color)]
end

function get_graph(grid_layout, system, cfg::CircleGraphCfg)
    ax = Axis(grid_layout[1, 1])
    circles =  [mesh!(ax, Circle(Point2f(0, 0), 1.0), color=:red, shading=NoShading)]
    drawn_borders(ax, system.space_cfg.geometry_cfg)
    CircleGraph(ax, circles)
end

function update_graph(graph::CircleGraph, system)
    data = get_graph_data(system, system.state)
    new_pos = Point2f[(data.pos[1, i], data.pos[2, i]) for i in 1:size(data.pos)[2]]
    update_circles!(graph, new_pos, data.radius, data.color)
end

end