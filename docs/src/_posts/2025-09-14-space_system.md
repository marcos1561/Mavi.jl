# Mavi's Space System
Here's a quick post explain how Mavi deals with systems in spaces with different geometries and border conditions.

# Defining a Space
A Space is defined by its geometry and how the borders (or walls) should behave, so there is a clear separation between geometry and wall behavior. Geometries are types derived from `GeometryCfg`, while walls are types derived from `WallType`. A space is simply the combination of a `GeometryCfg` and `WallType`, in code, this look like this

```julia
abstract type GeometryCfg end
abstract type WallType end

struct SpaceCfg{W<:WallType, G<:GeometryCfg} 
    wall_type::W
    geometry_cfg::G
end
```
when creating a system, one needs to provide some `SpaceCfg`.

For instance, we can create a rectangular geometry (with an arbitrary number of dimension) and periodic walls in the fallowing way

```julia
using StaticArrays # Needed to use SVector

# "Rectangle" with N dimensions using the number type T
struct RectangleCfg{N, T<:Number} <: GeometryCfg
    bottom_left::SVector{N, T}
    size::SVector{N, T} # size[i] is the rectangle's length in the i'th dimension
end

struct PeriodicWalls <: WallType end  
```

the actual code dealing with this kind of space should be in a specialized function, see next section.

# Dealing with Walls Behavior
Walls collisions should be resolved in the function `walls!(system, space_cfg::SpaceCfg)`, this function should dispatch on every kind of space there exists, for instance, with a rectangular geometry and periodic walls, the following method should be created.

```julia
function walls!(system, space_cfg::SpaceCfg{PeriodicWalls, RectangleCfg{N, T}}) where {N, T}
    # code to resolve walls collisions here
```

`walls!` must modify particles positions when necessary. Usually, `walls!` should be called after the update of positions by the equations of motion.

# Calculating Distances
Sometimes we need to know the information about the space to calculate distances (such as periodic boundaries with rectangle geometry), so there exists the function `calc_diff(r1, r2, space_cfg)`, which calculates the difference between positions `r1` and `r2`, and should dispatch on `space_cfg`. For instance,

```julia
# Default behavior
function calc_diff(r1, r2, space_cfg)
    @inbounds dr = r1 - r2
    return dr
end

# Periodic boundaries with rectangle geometry
function calc_diff(r1, r2, space_cfg::SpaceCfg{PeriodicWalls, G}) where G <: RectangleCfg
    size_vec = space_cfg.geometry_cfg.size
    dr = r1 - r2
    dr = dr - (abs.(dr) .> (size_vec / 2)) .* copysign.(size_vec, dr)
    return dr
end
```

# Using Multiple Spaces
It is possible to use multiple spaces, this is useful to insert obstacles in the simulation or create complex geometries using simple shapes. This functionality is accomplished by defining the types

```julia
struct ManyGeometries{G <: Tuple} <: GeometryCfg
    list::G
end

struct ManyWalls{W <: Tuple} <: WallType 
    list::W
end  
```

then, the walls behavior can simply be done like this

```julia
function walls!(system::System, space_cfg::SpaceCfg{W, G}) where {W <: ManyWalls, G <: ManyGeometries}
    for (wall_cfg, geom_cfg) in zip(space_cfg.wall_type.list, space_cfg.geometry_cfg.list)
        walls!(system, SpaceCfg(wall_cfg, geom_cfg))
    end
end
```

to make life easier for users, there exists a constructor for `SpaceCfg` to create composite spaces

```julia
# `spaces_cfg_list` must be a list of pairs of type (WallType, GeometryCfg).
function SpaceCfg(spaces_cfg_list)
    walls = []
    geometries = []

    for (w, g) in spaces_cfg_list
        push!(walls, w)
        push!(geometries, g)
    end

    SpaceCfg(
        wall_type=ManyWalls(Tuple(walls)),
        geometry_cfg=ManyGeometries(Tuple(geometries)),
    )
end
```

in fact, users should always use this constructor, because it ensures that the list of geometries and walls are in the correct order.
> OBS: Some functions in Mavi needs to know information about the space where particles are, therefore, when using multiple spaces, this information becomes not defied (because some spaces might just be obstacles). To mitigate this, the first space, in a space with multiple spaces, is considered as the "main" one.