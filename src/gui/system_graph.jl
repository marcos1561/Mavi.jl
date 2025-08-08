module SystemGraphs

export MainGraph, GraphCfg, GraphComp, GraphCompCfg, GraphCompDebug
export MainGraphCfg, CircleGraphCfg, ScatterGraphCfg
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

function drawn_borders(ax, space_cfg::RectangleCfg; adjust_lims=true)
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

    if adjust_lims
        dx = l*0.1
        dy = h*0.1
        xlims!(ax, x - dx, x + l + dx)
        ylims!(ax, y - dy, y + h + dy)
    end
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


function get_color_map(cmap::AbstractVector, num_types) 
    cmap_out = Vector{RGBf}(undef, num_types)
    for (idx, c) in enumerate(cmap)
        if c isa RGBf
            c_rgb = c
        elseif c == :random
            c_rgb = rand(RGBf)
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

function get_color_map(cmap::Symbol, num_types)
    try
        cmap = colorschemes[cmap]
    catch end

    if cmap == :random
        return [rand(RGBf) for _ in 1:num_types]
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

    # try
    #     cmap = colorschemes[cmap]
    # catch end

    # num_active = length(active_ids)

    # if cmap isa RGBf
    #     colors[1:num_active] .= cmap
    # elseif cmap isa Symbol 
    #     if cmap == :random 
    #         type_to_color = Dict(t => rand(RGBf) for t in unique(types))

    #         for (i, idx) in enumerate(active_ids)
    #             colors[i] = type_to_color[types[idx]]
    #         end
    #     else
    #         colors[1:num_active] .= RGBf(GLMakie.to_color(cmap))
    #     end
    # elseif cmap isa ColorScheme
    #     for idx in 1:num_active
    #         colors[idx] = cmap[rand(Float64)]
    #     end
    # else
    #     for (i, idx) in enumerate(active_ids)
    #         t = types[idx]
    #         if t > length(cmap) || t < 1
    #             colors[i] = RGBf(GLMakie.to_color(:pink))
    #         elseif cmap[t] == :random
    #             colors[i] = rand(RGBf)
    #         else
    #             colors[i] = cmap[t]
    #         end
    #     end
    # end
    # cv = @view colors[1:num_active]
    # return cv
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

function get_graph_data(graph_cfg::GraphCfg, system, state::State)
    return (pos=state.pos, types=Vector(1:system.num_p), num_p=system.num_p)
end
@inline get_graph_data(graph::GraphCfg, system) = get_graph_data(graph, system, system.state)

update_graph_data(graph::Graph, system, state::State) = DefaultDict{Symbol, Bool}(false)
@inline update_graph_data(graph::Graph, system) = update_graph_data(graph, system, system.state)


abstract type GraphCompCfg <: GraphCfg end
abstract type GraphComp <: Graph end
abstract type GraphCompDebug <: GraphComp end

get_graph_data(graph_cfg::GraphCompCfg, system, state::State) = collect(1:length(state.pos))

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
    
    types = get_graph_data(cfg, system, system.state)
    cmap = get_color_map(cfg.colors_map, length(types))
    colors = Vector{eltype(cmap)}(undef, num_total_particles)

    scatter_plot = scatter!(ax, [zero(eltype(pos_obs[]))]; cfg.kwargs...)

    comp = ScatterGraph(types, colors, cmap, scatter_plot, (), pos_obs, cfg)
    update_graph(comp, system)
    return comp

    # data = get_graph_data(cfg, system, system.state)
    # colors = Vector{RGBf}(undef, length(system.state.pos))

    # types_obs = Observable(types)
    # # colors_obs = lift(types_obs) do types
    # # colors_obs = lift(pos_obs) do _
    # #     # get_colors(cfg.colors_map, colors, types)
    # #     particles_ids = get_particles_ids(system.state, system.dynamic_cfg)
    # #     get_colors!(colors, types, cmap, particles_ids)
    # # end

    # pos_points_obs = lift(pos_obs) do pos
    #     [Point2f(p) for p in pos]
    # end
    
    # particles_ids = get_particles_ids(system.state, system.dynamic_cfg)
    # colors = get_colors!(colors, types, cmap, particles_ids)

    # scatter_plot = scatter!(ax, pos_points_obs; color=colors, cfg.kwargs...)
        
    # ScatterGraph(types, cmap, scatter_plot, (types_obs=types_obs,), pos_obs, cfg)
end

function update_graph(comp::ScatterGraph, system)
    update_list = get_comp_update_data(comp)(comp, system)
    notify_comp_observables(comp, update_list)

    particles_ids = get_particles_ids(system, system.state)
    colors = get_colors!(comp.colors, comp.types, comp.cmap, particles_ids)
    pos = comp.pos_obs[]
    points = [Point2f(p) for p in pos]
    Makie.update!(comp.plot, points; color=colors)
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

    types = get_graph_data(cfg, system, system.state)
    cmap = get_color_map(cfg.colors_map, length(types))
    colors = Vector{eltype(cmap)}(undef, num_total_particles)

    circles_plot = poly!(ax,
        [Circle(Point2f(0, 0), 1)],
        strokecolor=:black,
        strokewidth=0.5,
    )
    comp = CircleGraph(types, colors, radius, cmap, circles_plot, (), pos_obs, cfg)
    update_graph(comp, system)
    return comp
end

function update_graph(comp::CircleGraph, system)
    update_list = get_comp_update_data(comp)(comp, system)
    notify_comp_observables(comp, update_list)

    particles_ids = get_particles_ids(system, system.state)
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
        update_graph(c, system)
        # to_update = get_comp_update_data(c)(c, system)
        # update_all = to_update === nothing
        # for (name, obs) in pairs(get_comp_obs_list(c))
        #     if update_all || to_update[name]
        #         notify(obs)
        #     end
        # end
    end
    notify(graph.pos_obs)
end

end