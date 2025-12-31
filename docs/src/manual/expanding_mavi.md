# Expanding Mavi
Mavi was designed to be expanded. In this section, we will illustrate how to expand many of its parts.

## How to use custom forces?
Forces are calculated based on `DynamicCfg` with the function `calc_interaction(i, j, dynamic_cfg, system)` (which is in the module `Mavi.Integration` and in the file [`integration.jl`](https://github.com/marcos1561/Mavi.jl/blob/main/src/integration.jl)). So, one can use custom forces creating a new `DynamicCfg` and a method for `calc_interaction`

Example: Let's create a constant radial force (with modulus `force`) applied only when the distance between particles is less than `min_dist`

```julia
using Mavi.Configs
using Mavi.Integration

# New DynamicCfg struct
struct RadialForce <: DynamicCfg
    force::Float64
    min_dist::Float64
end

# Creating new method used in the forces calculation.
function Integration.calc_interaction(i, j, dynamic_cfg::RadialForce, system)
    pos = system.state.pos
    dr = calc_diff(pos[i], pos[j], system.space_cfg)
    dist = sqrt(sum(dr.^2))

    force = dynamic_cfg.force
    min_dist = dynamic_cfg.min_dist

    if dist > min_dist
        return zero(eltype(pos))
    end

    return force * dr / dist
end

# This method is used to draw the particles
Configs.particle_radius(dynamic_cfg::DynamicCfg) = dynamic_cfg.min_dist/2
```

After that, a `System` can be created as usual, with `dynamic_cfg` set to `RadialForce`. A complete example can be found here [`custom_force.jl`](https://github.com/marcos1561/Mavi.jl/blob/main/examples/custom_force.jl).

!!! note

    The function that computes all the pairwise forces is `calc_forces!`, it is inside the module `Mavi.Integration` and in the file [`integration.jl`](https://github.com/marcos1561/Mavi.jl/blob/main/src/integration.jl).

## How to use different equations of motion?
Just create your own step function with the appropriate equations of motion. Remember, you can reuse the general method to compute pairwise forces (`Mavi.Integration.calc_forces!`) and resolve wall collisions (`Mavi.Integration.walls!`), thus just focusing on the equations of motion.

Mavi already has functions to integrate some equations of motion, such as:
- `update_verlet!`: Newton Second Law integrated using the verlet method. 
- `update_rtp!`: Run-and-Tumble particles.
- `update_szabo!`: Szabo particles, which follow the equations of motion introduced in the paper ["Phase transition in the collective migration of tissue cells: experiment and model"](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.74.061908) by Szabó, B., et al. (Physical Review E, 2006).

if the state your equation of motion needs is not well described by the states Mavi already has, you can create a new state.

## How do wall collisions work?
The method `Mavi.Integration.walls!(system, space_cfg)` is responsible for resolving wall collisions. It dispatches on `space_cfg`, therefore, for instance, to implement rigid walls in a rectangular geometry (this method is already implemented on `Mavi.jl`), one should create the method

```julia
function walls!(system, space_cfg::SpaceCfg{RigidWalls, RectangleCfg})
    ...
end
```

`walls!` should directly change the state in `system`, resolving the collisions.

!!! note

    See more about how `Mavi` deals with wall collisions in the blog post [Mavi's Space System](@ref).
