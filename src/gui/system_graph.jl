module SystemGraphs

export Graph, GraphCfg, GraphComp, GraphCompCfg, GraphCompDebug
export update_graph, get_graph
export MainGraph
export 
    ManyGraphsCfg, MainGraphCfg, CircleGraphCfg, ScatterGraphCfg, 
    NumsGraphCfg, MovableObjectsCfg, TotalForceGraphCfg, ChunksGraphCfg
export drawn_borders, colors_from_cmap, get_graph_cfg

using GLMakie, ColorSchemes, DataStructures, Random, StaticArrays
using Mavi.Systems
using Mavi.States
using Mavi.Configs
using Mavi.MovableObjects
using Mavi.ChunksMod
using Mavi.Errors

"Drawn the borders of `geometry_cfg`."
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

# "Numbers of types given a list of "
# function get_num_types(types, cmap)
#     if typeof(cmap) <: AbstractArray
#         return max(length(types), length(cmap))
#     else
#         return length(types)
#     end
# end

"Returns `num` random colors from the colormap `cmap`."
function colors_from_cmap(cmap, num)
    if cmap isa Symbol
        cmap = colorschemes[cmap]
    end
    return Vector{Any}([cmap[rand(Float64)] for _ in 1:num])
end

"Convert `cmap` to a list of colors with `num_types` elements."
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
        return [cmap[rand(rng, Float32)] for _ in 1:num_types]            
    else
        color = RGBf(GLMakie.to_color(cmap))
        return fill(color, num_types)
    end
end

get_color_map(cmap::ColorScheme, num_types; rng=nothing) = cmap

"""
Returns a list of colors for the given particles, to be used in plotting.

# Arguments
- `colors`: 
    Preallocated vector to store colors for each particle.

- `types`: 
    Vector of indices for each particle.

- `cmap`: 
    Vector mapping types to colors.

- `particles_ids`: 
    Vector of particle indices.
"""
function get_colors!(colors, types, cmap, particles_ids)
    count = 1
    for idx in particles_ids
        colors[count] = cmap[types[idx]]
        count += 1
    end

    cs = @view colors[1:count-1]
    return cs
end

"""
Returns a list of radius for the given particles, to be used in plotting.

# Arguments
- `radius`: 
    Preallocated vector to store radius for each particle.

- `dynamic_cfg`: 
    System dynamic configuration.

- `state`: 
    System state.

- `particles_ids`: 
    Vector of particle indices.
"""
function get_radius!(radius, dynamic_cfg, state, particles_ids)
    count = 1
    for idx in particles_ids
        radius[count] = get_particle_radius(dynamic_cfg, state, idx)
        count += 1
    end
    rs = @view radius[1:count-1]
    return rs
end

# =
# Graphs
# = 

abstract type GraphCfg end
abstract type Graph end

"Construct a Graph for `ax` given its configurations `cfg`."
function get_graph(ax, system, cfg::GraphCfg) 
    throw(NotImplementedError(cfg))
    # @error "No graph constructor for $(typeof(cfg))"
end

"Updates the graph to the next frame."
function update_graph(graph::Graph, system) end


@inline get_graph_data(graph_cfg, system::System) = get_graph_data(graph_cfg, get_particles_state(system.state))

"""
Data used by Graph.

# Returns
    Named Tuple with two fields:

    - pos: Particles positions
    - types: Particles indices in the `pos` Vector.

    OBS: pos is a reference to system.state.pos
"""
function get_graph_data(graph_cfg::GraphCfg, state::State)
    num_p = length(state.pos)
    return (pos=state.pos, types=Vector(1:num_p))
end

@inline update_graph_data(graph, system::System) = update_graph_data(graph, get_particles_state(system.state))

"""
Updates the graph data but do not notify observables, only return a
dictionary of observable to be updated.

# Returns
There are two return options:

- DefaultDict{Symbol, Bool}  
    Dictionary mapping observable names (as symbols) to bools. True if
    the observable needs to be notified false otherwise.

- Nothing  
    Informs that every observable should be notified.
"""
update_graph_data(graph::Graph, state::State) = DefaultDict{Symbol, Bool}(false)

# =
# Graph Components
# = 

abstract type GraphCompCfg <: GraphCfg end
abstract type GraphComp <: Graph end
abstract type GraphCompDebug <: GraphComp end

"""
Data used by a component Graph.

# Returns
    Vector with particles indices in the pos Vector.
"""
get_graph_data(graph_cfg::GraphCompCfg, state::State) = Vector(1:length(state.pos))

"Returns the function to update the graph component data (particles types)."
get_comp_update_data(c::GraphComp) = c.cfg.update_data
get_comp_update_data(c::GraphCompDebug) = (_, _) -> DefaultDict{Symbol, Bool}(false)

"Returns the Named Tuple of observables in the graph component `c`"
get_comp_obs_list(c::GraphComp) = c.obs_list
get_comp_obs_list(c::GraphCompDebug) = ()

"Returns the default number of unique types used to color particles."
get_default_num_types(cfg::GraphCompCfg, state) = length(state.pos)

"Returns the number of unique types used to color particles."
function get_comp_num_types(cfg::GraphCompCfg, system)
    if typeof(cfg.colors_map) <: AbstractArray
        num_types = length(cfg.colors_map)
    else
        num_types = get_default_num_types(cfg, get_particles_state(system.state))
    end

    return num_types
end

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

get_graph_cfg(cfg) = cfg
get_graph_cfg(cfg::G) where {G<:GraphCompCfg} = MainGraphCfg(cfg)
get_graph_cfg(cfg::Union{Vector, Tuple}) = ManyGraphsCfg(Tuple(cfg))


struct ScatterGraphCfg{C, KW, F} <: GraphCompCfg
    colors_map::C
    kwargs::KW
    update_data::F
    rng::AbstractRNG
end

"""
Drawn particles using scatter.

# Arguments
- `colors_map`:  
    Mapping of types to colors. It can be:
    - :random  
        Random colors for each type
    - A colormap name as symbol  
        Random colors from the given colormap.
    - Vector of colors for each type.

- `kwargs`:  
    kwargs passed to scatter

- `update_data`:  
    Function to update the particles types. Defaults to no update.

- `rng`:  
    RNG object. Default to the global RNG.
"""
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

struct ScatterGraph{T, C, CMAP, P, O, PosObs} <: GraphComp
    types::Vector{T}
    colors::Vector{C}
    cmap::CMAP
    plot::P
    obs_list::O
    pos_obs::PosObs
    cfg::ScatterGraphCfg
end

function get_graph(ax, pos_obs, system, cfg::ScatterGraphCfg)
    num_total_particles = length(get_particles_state(system.state).pos)
    
    num_types = get_comp_num_types(cfg, system)

    types = get_graph_data(cfg, system)
    cmap = get_color_map(cfg.colors_map, num_types, rng=cfg.rng)
    colors = Vector{eltype(cmap)}(undef, num_total_particles)

    if cmap isa ColorScheme
        types = Float64.(types)
    end

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
    points = [Point(p...) for p in pos]

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

"""
Drawn particles as circles.

# Arguments
- `circle_radius`:  
    Radius used to draw the circles. If -1 (the default), infer particles radius 
    using system dynamic configurations. 

- `stroke_cfg`:  
    Stroke circle configurations.

- `colors_map`:  
    Mapping for types to colors. It can be:
    - :random  
        Random colors for each type
    - A colormap name as symbol  
        Random colors from the given colormap.
    - Vector of colors for each type.

- `update_data`:  
    Function to update the particles types. Defaults to no update.

- `rng`:  
    RNG object. Default to the global RNG.
"""
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
    num_total_particles = length(get_particles_state(system.state).pos)

    # TODO: Make cfg.circle_radius work with dispatch system
    if cfg.circle_radius == -1
        radius = Vector{Float64}(undef, num_total_particles)
    else
        if typeof(cfg.circle_radius) <: Number
            radius = fill(Float64(cfg.circle_radius), num_total_particles) 
        elseif length(cfg.circle_radius) != num_total_particles
            throw(ArgumentError(
                "Length of `cfg.circle_radius` ($(length(cfg.circle_radius))) should be " *
                "the same as number of particles ($(num_total_particles))"
            ))
        else
            radius = cfg.circle_radius
        end
    end

    num_types = get_comp_num_types(cfg, system)

    types = get_graph_data(cfg, system)
    cmap = get_color_map(cfg.colors_map, num_types, rng=cfg.rng)
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
    
    if comp.cfg.circle_radius == -1
        radius = get_radius!(comp.radius, system.dynamic_cfg, get_particles_state(system.state), particles_ids)
    else
        radius = comp.radius[1:length(particles_ids)]
    end

    pos = comp.pos_obs[]
    
    circles = [Circle(Point2f(p), r) for (p, r) in zip(pos, radius)]
    Makie.update!(comp.plot, circles; color=colors)
end

@kwdef struct NumsGraphCfg <: GraphCompCfg 
    offset=nothing
    kwargs = Dict()
end

struct NumsGrah{O, P} <: GraphComp 
    pos_obs::O
    plot::P
    cfg::NumsGraphCfg
end

function get_graph(ax, pos_obs, system, cfg::NumsGraphCfg)
    kwargs = cfg.kwargs

    if !(:align in keys(kwargs))
        kwargs[:align] = (:center, :center)
    end

    plot = text!(ax, [zero(eltype(pos_obs[]))]; cfg.kwargs...)
    graph = NumsGrah(pos_obs, plot, cfg)
    update_graph(graph, system)

    return graph
end

function update_graph(comp::NumsGrah, system)
    pos = comp.pos_obs[]
    particles_ids = get_particles_ids(system)
    
    offset = comp.cfg.offset
    if offset === nothing
        offset = zero(eltype(pos))
    end

    points = [Point2f(p + offset) for p in pos]
    text = [string(i) for i in particles_ids]

    Makie.update!(comp.plot, points; text=text)
end
    
"""
List of components to be drawn. Defaults to
drawn particles as circles.
"""
struct MainGraphCfg{T<:Tuple, G} <: GraphCfg
    comps_cfgs::T
    geometry_kwargs::G
end
function MainGraphCfg(comps=nothing; geometry_kwargs=nothing)
    if geometry_kwargs === nothing
        geometry_kwargs = ()
    end
    MainGraphCfg(get_main_graph_components(comps), geometry_kwargs)
end

get_main_graph_components(comps::Union{Vector, Tuple}) = tuple(comps...)
get_main_graph_components(comps::GraphCompCfg) = (comps,)
get_main_graph_components(comps) = get_main_graph_components(CircleGraphCfg())
    
struct MainGraph{A, P, O, C<:Tuple} <: Graph
    ax::A
    pos::P
    pos_obs::O
    components::C
    cfg::MainGraphCfg
end

"Build the Graph respective to `cfg` using the given `grid_layout` from Makie."
function get_graph(ax, system, cfg::MainGraphCfg)
    data = get_graph_data(cfg, system)
    pos_obs = Observable(data.pos)

    drawn_borders(ax, system.space_cfg.geometry_cfg; cfg.geometry_kwargs...)
    
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

function update_graph_data(graph::MainGraph, state::State)
    update_graph_data(graph, state.part_ids, state)
end

function update_graph(graph::MainGraph, system)
    update_graph_data(graph, system)
    for c in graph.components
        update_graph(c, system)
    end
    notify(graph.pos_obs)
end

# ==
# ManyGraphs
# ==

struct ManyGraphsCfg{T<:Tuple} <: GraphCfg
    graphs_cfg::T
end

struct ManyGraphs{T<:Tuple} <: Graph
    graphs::T
end

function get_graph(ax, system, cfg::ManyGraphsCfg)
    ManyGraphs(tuple([get_graph(ax, system, get_graph_cfg(cfg_i)) for cfg_i in cfg.graphs_cfg]...))
end

function update_graph(graph::ManyGraphs, system)
    for graph_i in graph.graphs
        update_graph(graph_i, system)
    end
end

# ==
# MovableObjects
# ==
struct MovableObjectsCfg <: GraphCfg end

struct MovableObjectsGraph{M <: MovableObject} <: Graph 
    cfg::MovableObjectsCfg
    movable_objects::Vector{M}
    objects_states::Vector
end

function init_movable_object_graph(state::LineState, object::MovableObject, ax)
    line_obs = Observable(state)
    lines!(ax, lift(l -> [l.p1, l.p2], line_obs), color="black")

    return line_obs
end

function update_movable_object_graph(graph_state::Observable, object)
    graph_state[] = object.state
end

function get_graph(ax, system, cfg::MovableObjectsCfg)
    graph_states = []
    movable_objects = get_movable_objects(system.state)
    for object in movable_objects
        push!(graph_states, init_movable_object_graph(object.state, object, ax)) 
    end
    MovableObjectsGraph(cfg, movable_objects, graph_states)
end

function update_graph(graph::MovableObjectsGraph, system)
    for (graph_state, object) in zip(graph.objects_states, graph.movable_objects)
        update_movable_object_graph(graph_state, object)
    end
end

# = 
# Forces
# =
@kwdef struct TotalForceGraphCfg{C} <: GraphCompCfg
    color::C = "black"
    lengthscale::Float64 = 1.0
end

struct TotalForceGraph{P, C} <: GraphCompDebug
    arrows::P
    cfg::C
end

function SystemGraphs.get_graph(ax, pos_obs, system, cfg::TotalForceGraphCfg)
    plot = arrows2d!(ax, [0.0], [0.0], [0.0], [0.0], lengthscale=cfg.lengthscale)
    graph = TotalForceGraph(plot, cfg)
    SystemGraphs.update_graph(graph, system)
    return graph
end

function SystemGraphs.update_graph(graph::TotalForceGraph, system)
    pos = system.state.pos
    p_forces = get_forces(system)
    
    n = length(get_particles_ids(system))
    x = Vector{Float64}(undef, n)
    y = Vector{Float64}(undef, n)
    u = Vector{Float64}(undef, n)
    v = Vector{Float64}(undef, n)
    
    idx = 1
    for i in get_particles_ids(system)
        x[idx] = pos[i][1]
        y[idx] = pos[i][2]
        u[idx] = p_forces[i][1]
        v[idx] = p_forces[i][2]
        idx += 1
    end

    Makie.update!(graph.arrows, x, y, u, v)
end

# ==
# Chunks
# ==

@kwdef struct ChunksGraphCfg{C} <: GraphCfg 
    color::C="gray"
    see_particle_chunk=true
end

struct ChunksGraph{P} <: Graph 
    cfg::ChunksGraphCfg
    p_chunk_plot::P
end

function get_graph(ax, system, cfg::ChunksGraphCfg)
    chunks = system.chunks
    if chunks === nothing
        return nothing
    end

    bbox = chunks.space_cfg.geometry_cfg
    bl = bbox.bottom_left
    l, h = chunks.chunk_length, chunks.chunk_height
    for i in 1:(chunks.num_cols + 1)
        x = bl.x + l * (i - 1)
        lines!(ax, [x, x], [bl.y, bl.y + bbox.height], color=cfg.color)
    end
    for i in 1:(chunks.num_rows + 1)
        y = bl.y + h * (i - 1)
        lines!(ax, [bl.x, bl.x + bbox.length], [y, y], color=cfg.color)
    end

    p_chunk_plot = nothing
    if cfg.see_particle_chunk
        p_chunk_plot = text!(ax, [zero(eltype(system.state.pos))])
    end
    graph = ChunksGraph(cfg, p_chunk_plot)
    update_graph(graph, system)
    graph
end

function update_graph(graph::ChunksGraph, system) 
    if graph.p_chunk_plot === nothing
        return
    end

    pos = system.state.pos
    chunks = system.chunks
    nc = chunks.num_cols
    nr = chunks.num_rows

    points = Point2f[]
    text = String[]
    for i in 1:nr, j in 1:nc
        for pid in get_chunk_particles(chunks, i, j)
            push!(points, pos[pid])
            push!(text, "($i, $j)")
        end
    end

    Makie.update!(graph.p_chunk_plot, points; text=text)
end

include("../rings/view.jl")

end