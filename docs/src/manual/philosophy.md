# Mavi Philosophy 
Every system of particles is represented by the same struct [`System`](https://github.com/marcos1561/Mavi.jl/blob/main/src/systems.jl), whose main members are:

- **state**: Data representing the particles' states.
- **space_cfg**: Information about the geometry of the space (including boundary behavior) where the particles are located.
- **dynamic_cfg**: Configuration of the dynamics between particles.
- **int_cfg**: Configuration of how to integrate the equations of motion.

In order to make things work, all these members' types are parametric, so Julia's multiple dispatch mechanism can shine.

A `System` should be integrated by a function that executes one time step at a time. Usually, a step function looks like this:

```julia
function step!(system)
    update!(system)       # Integrates the equations of motion
    walls!(system)        # Resolves wall collisions
    update_time!(system)  # Updates time
    clean_forces!(system) # Sets forces array to zero
    calc_forces!(system)  # Calculates all the forces
end
```

!!! note

    Note that one needs to setup the system (calculate the forces and everything else needed in order to update the state) before start calling the `step!` function. This is the recommended approach, because after a `step!` execution, the system is in a state that reflects the current time, otherwise the system would be one time step behind.

Mavi already has some step functions (see [`integration.jl`](https://github.com/marcos1561/Mavi.jl/blob/main/src/integration.jl)), but users are free to create their own.

In the preceding sections we will talk about the main system's members in more detail: what they represent and how to create them.  

# State
The state of a system is defined by the equations that govern its dynamics. For example, an ordinary, linear, second-order differential equation in $x(t)$ requires that two variables be stored in memory (the variable $x$ itself and its derivative $\dot x$).
States must be subtypes of `State{N, T}` (see [`states.jl`](https://github.com/marcos1561/Mavi.jl/blob/main/src/states.jl) for examples), where `N` is the spacial dimension and `T` the numeric type (`Int32`, `Float64`, etc). Currently, two states are implemented in the core `Mavi.States` module:

- `SecondLawState`: has `pos` and `vel` as members, and both are of type `Vector{SVector}`. This state is useful when the equation to be solved is Newton's second law.

- `SelfPropelledState`: a state related to a system of particles with overdamped dynamics and a self-alignment mechanism.

The following example creates a state with 3 particles using `SecondLawState`.

```julia
using Mavi.States

state = SecondLawState(
    pos = [1 2 3; 1 2 3],
    vel = [1.1 0.4 3; 1.2 2 -0.4],
)
```
!!! note

    Although `pos` and `vel` are stored as `Vector{SVector}`, they can be constructed with other types, such as matrices, as shown in the exemple above. 

## SpaceCfg
The space of a system is described by its geometric shape and how its boundaries behave:

- **Space geometry**: These are structures that inherit from `GeometryCfg` and contain all the necessary information to describe the geometry of the space. Examples: `RectangleCfg` and `CircleCfg`. 

- **Boundary behavior**: These are structures that inherit from `WallType` and do not necessarily need to have any fields; they simply indicate the boundary type. Examples: `RigidWalls` and `PeriodicWalls`.

`SpaceCfg` encapsulates these two structures. See [`configs.jl`](https://github.com/marcos1561/Mavi.jl/blob/main/src/configs.jl) for more details. It is possible to combine many `SpaceCfg` to form more complex spaces, this is done passing a list of `Tuple{WallType, GeometryCfg}` to  the `SpaceCfg` constructor (in this case, the first `SpaceCfg` is considered as the "main" one, and the others can be view as obstacles or holes).

Example: Configuring a rectangle with its bottom-left corner at the origin (default behavior) and with periodic boundaries:

```julia
using Mavi.Configs

space_cfg = SpaceCfg(
    geometry_cfg=RectangleCfg(
        length=5,
        height=10,
    ),
    wall_type=PeriodicWalls(),
)
```

## DynamicCfg
This object should have configurations about the interactions between particles, such as potential parameters. Every instance is a subtype of `DynamicCfg`.

As examples, the `Mavi.Configs` module has a special kind of `DynamicCfg` that represents potential forces, the `PotentialCfg`, examples are:

- Harmonic Truncated Potential: `HarmTruncCfg`
- Lennard-Jones Potential: `LenJonesCfg`

See more code details in [configs.jl](https://github.com/marcos1561/Mavi.jl/blob/main/src/configs.jl).

`DynamicCfg` is used to calculate forces, users should only implement a method for `Mavi.Integration.calc_interaction(i, j, dynamic_cfg, system)`, where i and j are particle indices (more info in its documentation [`integration.jl`](https://github.com/marcos1561/Mavi.jl/blob/main/src/integration.jl)). If using a `PotencialCfg`, users only need to implement a method for `Mavi.Configs.potential_force(dr, dist, potential::PotencialCfg)` that returns the force in the position `dr`.

Also, particles radius should be inferred from these configurations, so there exists a function `particle_radius(dynamic_cfg)`. 

## IntCfg
Configurations for how the equations of motion should be integrated. Every instance is a subtype of `AbstractIntCfg`.
The default integration configuration inside `Mavi.Configs` is `IntCfg`, which has two main members:
- **dt**: Time step used in the integration.
- **chunks_cfg**: Configurations of the technique used to speed distance calculations: the system is divided into boxes, or chunks, and the interaction between particles is calculated only for neighboring boxes, greatly reducing the simulation time for short range interactions. If its value is `nothing`, then this technique is not used.

The following example creates a system with a circular space, using the Lennard-Jones potential and simple time integration:

```julia
using Mavi
using Mavi.Configs

system = System(
    state=state, # assuming the state has already been created
    space_cfg=SpaceCfg(
        geometry_cfg=CircleCfg(radius=10, center=[0, 0]),
        wall_type=RigidWalls(),
    ), 
    dynamic_cfg=LenJonesCfg(sigma=2,epsilon=4),
    int_cfg=IntCfg(dt=0.01),
)
```

# About the step function
The step function should evolve a `System` in one time step. The actions usually needed to perform a time step are: 
- integrate equations of motion
- resolve wall collisions    
- calculate forces

Mavi already has a general system to compute pairwise interactions, some ready to use functions to integrate equations of motion (such as Newton's Second Law) and wall collisions resolvers, but users are free to create their own methods as specified in the next sections of this documentation.