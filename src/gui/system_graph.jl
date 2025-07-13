module SystemGraphs

export GraphCfg, DefaultGraphCfg
export CircleGraph, CircleGraphCfg

using GLMakie
using Mavi.States: State
using Mavi.Systems: System
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
    x, y = space_cfg.bottom_left
    lines!(ax, [x, y], [x+l, y], color=:black)
    lines!(ax, [x+l, y], [x+l, y+h], color=:black)
    lines!(ax, [x+l, y+h], [x, y+h], color=:black)
    lines!(ax, [y+h, x], [x, y], color=:black)

    dx = l*0.1
    dy = h*0.1
    xlims!(ax, x - dx, x + l + dx)
    ylims!(ax, y - dy, y + h + dy)
end

function drawn_borders(ax, space_cfg::LinesCfg)
    for line in space_cfg.lines
        p1, p2 = line.p1, line.p2
        lines!(ax, [p1.x, p2.x], [p1.y, p2.y], color=:black)
        # l, h = space_cfg.length, space_cfg.height
        # x, y = space_cfg.bottom_left
        # lines!(ax, [x, y], [x+l, y], color=:black)
        # lines!(ax, [x+l, y], [x+l, y+h], color=:black)
        # lines!(ax, [x+l, y+h], [x, y+h], color=:black)
        # lines!(ax, [y+h, x], [x, y], color=:black)

        # dx = l*0.1
        # dy = h*0.1
        # xlims!(ax, x - dx, x + l + dx)
        # ylims!(ax, y - dy, y + h + dy)
    end
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
    colors_map = [RGBf(GLMakie.to_color(:black))]
end
function DefaultGraphCfg(circle_radius, circle_rel, color_map::Vector{Symbol})
    color_map = RGBf.(GLMakie.to_color.(color_map))
    DefaultGraphCfg(circle_radius, circle_rel, color_map)
end

struct DefaultGraph{T1, T2, T3}
    ax::T1
    pos_obs::T2
    colors_obs::T3
    cfg::DefaultGraphCfg
end

function get_graph_data(system, state) end

function get_graph_data(graph_cfg::DefaultGraphCfg, system, state)
    return (pos=state.pos, types=fill(1, size(state.pos)[2]))
end

function get_colors(graph_cfg::DefaultGraphCfg, colors, types)
    count = 1
    for t in types
        if graph_cfg.colors_map[t] == :random
            colors[count] = rand(RGBf)
        else
            colors[count] = graph_cfg.colors_map[t]
        end
        count += 1
    end

    cv = @view colors[1:length(types)]
    return cv
end

"Build the Graph respective to `cfg` using the given `grid_layout` from Makie."
function get_graph(grid_layout, system, cfg::DefaultGraphCfg)
    ax = Axis(grid_layout[1, 1], aspect=DataAspect())

    data = get_graph_data(cfg, system, system.state)

    pos_obs = Observable(data.pos)

    colors = Vector{RGBf}(undef, size(data.pos, 2))
    colors_obs = Observable(get_colors(cfg, colors, data.types))

    drawn_borders(ax, system.space_cfg.geometry_cfg)
    
    if cfg.circle_radius == -1
        radius = particle_radius(system.dynamic_cfg) 
    else
        radius = cfg.circle_radius 
    end
    
    function points_to_circle(pos)
        [Circle(Point2f(pos[1, i], pos[2, i]), radius) for i in axes(pos, 2)]
    end
    
    # scatter!(ax, x_obs, y_obs, color=colors_obs)
    
    # Circles
    poly!(ax,
        lift(pos_obs) do pos
            points_to_circle(pos)
        end,
        # color = :dodgerblue,
        color = colors_obs,
        strokecolor = :black,
        strokewidth = 0.5,
    )

    DefaultGraph(ax, pos_obs, colors_obs, cfg)
end

function update_graph(graph::DefaultGraph, system)
    get_graph_data(graph.cfg, system, system.state)
    notify(graph.pos_obs)
end

"""
Drawn every particle as circle with given radius and color.
"""
struct CircleGraphCfg <: GraphCfg end

mutable struct CircleGraph{T1, T2}
    ax::T1
    circles::T2
end

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