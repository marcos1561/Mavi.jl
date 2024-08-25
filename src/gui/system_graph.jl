module SystemGraphs

export GraphCfg, DefaultGraphCfg

using GLMakie
using Mavi: State
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
end

function drawn_borders(ax, space_cfg::CircleCfg)
    arc!(ax, Point2f(0), space_cfg.radius, -π, π, color=:black)
end

function update_circles!(circle_points, circle_base, state::State)
    num_edges = length(circle_base)
    num_p = size(state.pos)[2]
    for i in 1:num_p
        begin_id = num_edges * (i-1) + 1
        end_id = begin_id + num_edges - 1 
        center = [state.pos[1, i], state.pos[2, i]]
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
end

struct DefaultGraph{T1, T2}
    ax::T1
    circle_base::Vector{Point2f}
    circle_points::Vector{Point2f}
    circle_points_obs::Observable{Vector{Point2f}}
    pos_obs::T2
end

"Build the Graph respective to `cfg` using the given `grid_layout` from Makie."
function get_graph(grid_layout, system, cfg::DefaultGraphCfg)
    ax = Axis(grid_layout[1, 1])

    x_obs = Observable(@view system.state.pos[1, :])
    y_obs = Observable(@view system.state.pos[2, :])
    pos_obs = [x_obs, y_obs]

    scatter!(ax, pos_obs...)
    # ax.aspect = 1

    drawn_borders(ax, system.space_cfg)

    if cfg.circle_radius == -1
        radius = particle_radius(system.dynamic_cfg) 
    else
        radius = cfg.circle_radius 
    end

    # Circles
    circle_base = get_circle(radius, cfg.circle_rel)
    circle_points = Vector{Point2f}(undef, system.num_p * length(circle_base))

    update_circles!(circle_points, circle_base, system.state)
    circle_points_obs = Observable(circle_points)
    linesegments!(ax, circle_points_obs, color=:black)

    DefaultGraph(ax, circle_base, circle_points, circle_points_obs, pos_obs)
end

function update_graph(graph::DefaultGraph, state)
    update_circles!(graph.circle_points, graph.circle_base, state)
    notify(graph.circle_points_obs)

    notify(graph.pos_obs[1])
    notify(graph.pos_obs[2])
end

end