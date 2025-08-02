module SystemGraphs

export GraphCfg, MainGraphCfg, CircleCfg, ScatterGraphCfg
export MainGraph, CircleGraphCfg, ScatterGraphCfg
export drawn_borders

using GLMakie, ColorSchemes, DataStructures
using Mavi.Systems
using Mavi.States: State
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
    # lines!(ax, [x, y], [x+l, y], color=:black)
    # lines!(ax, [x+l, y], [x+l, y+h], color=:black)
    # lines!(ax, [x+l, y+h], [x, y+h], color=:black)
    # lines!(ax, [y+h, x], [x, y], color=:black)
    
    lines!(ax, [x, x+l], [y, y], color=:black)
    lines!(ax, [x+l, x+l], [y+h, y], color=:black)
    lines!(ax, [x+l, x], [y+h, y+h], color=:black)
    lines!(ax, [x, x], [y+h, y], color=:black)

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
abstract type Graph end
abstract type GraphCompCfg <: GraphCfg end
abstract type GraphComp <: Graph end

struct ScatterGraphCfg{C, KW, F} <: GraphCompCfg
    colors_map::C
    kwargs::KW
    update_data::F
end
function ScatterGraphCfg(;colors_map=:random, update_data=nothing, kwargs=nothing)
    if colors_map isa Vector{Symbol}
        colors_map = RGBf.(GLMakie.to_color.(colors_map))
    end

    if update_data === nothing
        update_data = update_graph_data
    end

    if kwargs === nothing
        kwargs = ()
    end

    ScatterGraphCfg(colors_map, kwargs, update_data)
end

struct ScatterGraph{O} <: GraphComp
    types::Vector{Int}
    obs_list::O
    cfg::ScatterGraphCfg
end

function get_graph(ax, pos_obs, system, cfg::ScatterGraphCfg)
    data = get_graph_data(cfg, system, system.state)
    colors = Vector{RGBf}(undef, size(system.state.pos, 2))

    types_obs = Observable(data.types)
    colors_obs = lift(types_obs) do types
        get_colors(cfg.colors_map, colors, types)
    end
    
    scatter!(ax, pos_obs; color=colors_obs, cfg.kwargs...)
        
    ScatterGraph(data.types, (types_obs=types_obs,), cfg)
end

struct CircleGraphCfg{C, F} <: GraphCompCfg
    circle_radius::Float64
    circle_rel::Int
    colors_map::C
    update_data::F
end
function CircleGraphCfg(;circle_radius=-1.0, circle_rel=20, colors_map=:random, 
    update_data=nothing)
    # if colors_map === nothing
    #     colors_map = RGBf(GLMakie.to_color(:black))
    # end

    if colors_map isa Vector{Symbol}
        colors_map = RGBf.(GLMakie.to_color.(colors_map))
    end

    if update_data === nothing
        update_data = update_graph_data
    end

    CircleGraphCfg(circle_radius, circle_rel, colors_map, update_data)
end

struct CircleGraph{C, O} <: GraphComp
    types::Vector{Int}
    colors::Vector{C}
    obs_list::O
    cfg::CircleGraphCfg
end

function get_graph(ax, pos_obs, system, cfg::CircleGraphCfg)
    if cfg.circle_radius == -1
        radius = particles_radius(system, system.dynamic_cfg)
    else
        if cfg.circle_radius <: Number
            radius = fill(cfg.circle_radius, size(data.pos, 2)) 
        elseif length(cfg.circle_radius) != size(data.pos, 2)
            throw(ArgumentError(
                "Length of `cfg.circle_radius` ($(length(cfg.circle_radius))) should be " *
                "the same as number of particles ($(size(data.pos, 2)))"
            ))
        end
    end

    data = get_graph_data(cfg, system, system.state)
    colors = Vector{RGBf}(undef, size(system.state.pos, 2))
    # colors = Vector{RGBf}(undef, length(system.state.pos))

    types_obs = Observable(data.types)
    colors_obs = lift(types_obs) do types
        get_colors(cfg.colors_map, colors, types[1:get_num_total_particles(system)])
    end

    function points_to_circle(pos)
        # [Circle(Point2f(pos[1, i], pos[2, i]), radius[i]) for i in axes(pos, 2)]
        [Circle(Point2f(pos[i]), radius[i]) for i in eachindex(pos)]
    end

    poly!(ax,
        lift(pos_obs) do pos
            points_to_circle(pos)
        end,
        # color = :dodgerblue,
        color=colors_obs,
        strokecolor=:black,
        strokewidth=0.5,
    )

    CircleGraph(data.types, colors, (types_obs=types_obs,), cfg)
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

function get_colors(cmap, colors, types)
    try
        cmap = colorschemes[cmap]
    catch end

    if cmap isa RGBf
        colors .= cmap
    elseif cmap isa Symbol 
        if cmap == :random 
            type_to_color = Dict(t => rand(RGBf) for t in unique(types))

            for idx in eachindex(colors)
                colors[idx] = type_to_color[types[idx]]
            end
        else
            colors .= RGBf(GLMakie.to_color(cmap))
        end
    elseif cmap isa ColorScheme
        for idx in eachindex(colors)
            colors[idx] = cmap[rand(Float64)]
        end
    else
        count = 1
        for t in types
            if t > length(cmap) || t < 1
                colors[count] = RGBf(GLMakie.to_color(:pink))
            elseif cmap[t] == :random
                colors[count] = rand(RGBf)
            else
                colors[count] = cmap[t]
            end
            count += 1
        end
    end
    cv = @view colors[1:length(types)]
    return cv
end

function get_graph_data(graph_cfg::GraphCfg, system, state::State)
    return (pos=state.pos, types=Vector(1:system.num_p), num_p=system.num_p)
end

update_graph_data(graph::Graph, system, state::State) = DefaultDict{Symbol, Bool}(false)

@inline get_graph_data(graph::Graph, system) = get_graph_data(graph, system, system.state)
@inline update_graph_data(graph::Graph, system) = update_graph_data(graph, system, system.state)

"Build the Graph respective to `cfg` using the given `grid_layout` from Makie."
function get_graph(ax, system, cfg::MainGraphCfg)
    data = get_graph_data(cfg, system, system.state)    
    # pos_view = @view data.pos[:, 1:data.num_p]
    pos_view = @view data.pos[1:data.num_p]
    pos_obs = Observable(pos_view)

    drawn_borders(ax, system.space_cfg.geometry_cfg)
    
    graphs = Vector{GraphComp}()
    for g in cfg.comps_cfgs
        push!(graphs, get_graph(ax, pos_obs, system, g))
    end

    # scatter!(ax, x_obs, y_obs, color=colors_obs)
    
    # Circles
    # poly!(ax,
    #     lift(pos_obs) do pos
    #         points_to_circle(pos)
    #     end,
    #     # color = :dodgerblue,
    #     color = colors_obs,
    #     strokecolor = :black,
    #     strokewidth = 0.5,
    # )

    # scatter_graph = nothing
    # if !(scatter_cfg === nothing)
    #     colors = Vector{RGBf}(undef, size(data.pos, 2))
    #     sc_colors_obs = Observable(get_colors(cfg.scatter_cfg.colors_map, colors, data.types))
    #     scatter!(pos_obs, color=colors_obs)
        
    #     scatter_graph = ScatterGraph(sc_colors_obs, cfg.scatter_cfg)
    # end

    graph = MainGraph(ax, data.pos, pos_obs, tuple(graphs...), cfg)
    update_graph(graph, system)
    graph
end

function update_graph(graph::MainGraph, system)
    update_graph_data(graph, system)
    for c in graph.components
        to_update = c.cfg.update_data(c, system)
        update_all = to_update === nothing
        for (name, obs) in pairs(c.obs_list)
            if update_all || to_update[name]
                notify(obs)
            end
        end
    end
    notify(graph.pos_obs)
end

end