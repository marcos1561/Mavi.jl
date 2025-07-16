# Mavi.jl
Mavi is a _Particle Dynamics Engine_.

Its goal is to provide a common structure for implementing particle dynamics simulations, allowing users to use default behaviors or create their own as needed.

Here are some `Mavi.jl` features:

- Pairwise interactions, such as Lennard-Jones or Harmonic-Truncated potentials.

- Spaces with different geometries, such as rectangular or circular.

- Spaces with different wall types, such as rigid or periodic walls.

- Quantity calculators, such as kinetic and potential energy.

- Visualization:

    - Real-time rendering of simulations using [Makie](https://docs.makie.org/v0.21/).

    - Video generation of simulations.

In addition to the core `Mavi` module, there is a module named `Mavi.Rings`, which implements active rings — a model used in this paper: ["Segregation in Binary Mixture with Differential Contraction among Active Rings"](https://link.aps.org/doi/10.1103/PhysRevLett.134.138401) by Teixeira, E., et al (Physical Review Letters, 2025).

# Installation
The easiest way to install Mavi is:

* Open the Julia REPL and enter package manager mode (by pressing `]` at the empty REPL prompt)
* Execute the command:
    ```julia
    pkg> dev https://github.com/marcos1561/Mavi.jl.git 
    ```
This will install Mavi in the currently active environment, placing its files in `~/.julia/dev/Mavi`.

Note: If you just want to use Mavi (not modify it), use `add` instead of `dev`.

# Mavi Philosophy 
Every system of particles is represented by the same struct [`System`](src/systems.jl), whose main members are:

- state: Data representing the particles' states.
- space_cfg: Information about the geometry of the space (including boundary behavior) where the particles are located.
- dynamic_cfg: Configuration of the dynamics between particles.
- int_cfg: Configuration of how to integrate the equations of motion.

In order to make things work, all these members' types are parametric, so Julia's multiple dispatch mechanism can shine.

A `System` should be integrated by a function that executes one time step at a time. Usually, a step function looks like this:

```julia
function step!(system)
    calc_forces!(system) # Calculates all the forces
    update!(system)      # Integrates equations of motion
    walls!(system)       # Resolves wall collisions
end
```

Mavi already has some step functions (see [`integration.jl`](src/integration.jl)), but users are free to create their own.

# About System Members
In this section we will talk about the main system's members: what they represent and how to create them.  

## State
The state of a system is defined by the equations that govern its dynamics. For example, an ordinary, linear, second-order differential equation in $x(t)$ requires that two variables be stored in memory (the variable $x$ itself and its derivative $\dot x$).
States must be subtypes of `State{T}` (see [`states.jl`](src/states.jl) for examples), where T is the type of the numbers representing the state (Int32, Float64, etc). Currently, two states are implemented in the core `Mavi.States` module:

`SecondLawState`: has pos and vel as members, both of which must be $2 \times N$ matrices. This state is useful when the equation to be solved is Newton's second law.
`SelfPropelledState`: a state related to a system of particles with overdamped dynamics and a self-alignment mechanism.

The following example creates a state with 3 particles using `SecondLawState`.

```julia
using Mavi.States

state = SecondLawState(
    pos = [1 2 3; 1 2 3],
    vel = [1.1 0.4 3; 1.2 2 -0.4],
)
```

## SpaceCfg
The space of a system is described by its geometric shape and how its boundaries behave:

- Space geometry: These are structures that inherit from `GeometryCfg` and contain all the necessary information to describe the geometry of the space. Examples: `RectangleCfg` and `CircleCfg`.

Boundary behavior: These are structures that inherit from `WallType` and do not necessarily need to have any fields; they simply indicate the type of boundary. Examples: `RigidWalls` and `PeriodicWalls`.

`SpaceCfg` encapsulates these two structures. See [`configs.jl`](src/configs.jl) for more details.

Example: Configuring a rectangle with its bottom-left corner at the origin and with periodic boundaries

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

As examples, `Mavi.Configs` module has these potential structures:

- Harmonic Truncated Potential: `HarmTruncCfg(ko, ro, ra)`
- Lennard-Jones Potential: `LenJonesCfg(sigma, epsilon)`

Details about parameters can be found in [configs.jl](src/configs.jl).

`DynamicCfg` is used to calculate forces, users should only implement a method for `calc_interaction(i, j, state, dynamic_cfg, space_cfg)` (inside module `Mavi.Integration`), where i and j are particles indices (more info in its documentation [`integration.jl`](src/integration.jl)).

Also, particles radius should be inferred from these configurations, so there exists a function `particle_radius(dynamic_cfg)`. 

## IntCfg
Configurations for how the equations of motion should be integrated. Every instance is a subtype of `AbstractIntCfg`.
The default integration configuration inside `Mavi.Configs` is `IntCfg`, which has two members:
- dt: Time step used in the integration.
- chunks_cfg: Configurations of the technique used to speed distance calculations: the system is divided into boxes, or chunks, and the interaction between particles is calculated only for neighboring boxes, greatly reducing the simulation time for short range interactions. If its value is `nothing`, then this technique is not used.

The following example creates a system with a circular space, using the Lennard-Jones potential and simple time integration:

```julia
using Mavi
using Mavi.Configs

system = System(
    state=state, # assuming the state has already been created
    space_cfg=CircleCfg(radius=10), 
    dynamic_cfg=LenJonesCfg(sigma=2,epsilon=4),
    int_cfg=IntCfg(dt=0.01),
)
```

# About the step function
The step function should evolve a `System` in one time step. The actions usually needed to perform a time step are: 
- calculate forces
- integrate equations of motion
- resolve wall collisions    

Mavi already has a general system to compute pair interactions, some ready to use functions to integrate equations of motion (such as Newton's Second Law) and wall collisions resolvers.

## How to use custom forces?
Forces are calculated based on `DynamicCfg` with the function `calc_interaction(i, j, state, dynamic_cfg, space_cfg)` (which is in the module `Mavi.Integration` and in the file [`integration.jl`](src/integration.jl)). So, one can use custom forces creating a new `DynamicCfg` and a method for `calc_interaction`

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
function Integration.calc_interaction(i, j, state, dynamic_cfg::RadialForce, space_cfg)
    dx, dy, dist = calc_diff_and_dist(i, j, state.pos, space_cfg)

    force = dynamic_cfg.force
    min_dist = dynamic_cfg.min_dist

    if dist > min_dist
        return 0.0, 0.0
    end

    # Force components
    fx = force * dx / dist
    fy = force * dy / dist

    return fx, fy
end

# This method is used to draw the particles
Configs.particle_radius(dynamic_cfg::DynamicCfg) = dynamic_cfg.min_dist/2
```

After that, a `System` can be created as usual, with `dynamic_cfg` set to `RadialForce`. A complete example can be found here [`custom_force.jl`](examples/custom_force.jl).

>Note: The function that computes all the pairwise forces is `calc_forces!()`, it is inside the module `Mavi.Integration` and in the file [`integration.jl`](src/integration.jl).


## How to use different equations of motion?
Just create your own step function with the appropriate equations' of motion. Remember, you can reuse the general method to compute pairwise forces (`calc_forces!()`) and wall collision detections, thus just focusing on the equations of motion.

Mavi already has functions to integrate some equations of motion, such as:
- `update_verlet!`: Newton Second Law integrated using the verlet method. 
- `update_rtp!`: Run-and-Tumble particles.
- `update_szabo!`: Szabo particles, which follow the equations of motion introduced in the paper ["Phase transition in the collective migration of tissue cells: experiment and model"](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.74.061908) by Szabó, B., et al. (Physical Review E, 2006).

## How do wall collisions work?
The method `walls!(system, space_cfg)` (inside module `Mavi.Integration`) is responsible to resolve wall collisions. It dispatches on `space_cfg`, therefore, for example, to implement rigid walls in a retangular geometry (this method is already implemented on `Mavi.jl`), one should create the method

```julia
function walls!(system, space_cfg::SpaceCfg{RigidWalls, RectangleCfg})
    ...
end
```

`walls!` should directly change the state in `system`, resolving the collisions.

# Physical Quantities
The `Quantities` module contains functions that return physical quantities of interest. Currently, two of them are defined:

- **Kinetic energy:**  
    The system's kinetic energy is calculated by the function `kinetic_energy(state)`, which takes the system's state as a parameter. Since it depends only on the velocities, it is independent of the chosen dynamic configuration.

- **Potential energy:**  
    The system's potential energy is calculated by the function `potential_energy(system, dynamic_cfg)`, which takes the system and the dynamic configuration as arguments, since it depends on the chosen interaction between particles.  
    The user should ensure that the correct dispatch for the `potential_energy()` function is defined for any custom dynamic configuration.

The example [print_energy.jl](examples/print_energy.jl) uses these functions.

# Visual Interface
Mavi has a visual interface (built entirely with [Makie](https://docs.makie.org/v0.21/)) whose purpose is to serve as a visual debugging tool for the system being explored. The UI structure essentially has two main elements:

1. A plot where the system's particles are rendered in real time.
2. An information panel showing the state of the system and the execution of the program.

<!-- ![alt text](docs/images/ui_components.png "Title") -->
<img src="docs/images/ui_components.png" alt="Componentes da UI" width="500"/>

## Animating the system
Given that a `step!(system)` function has already been implemented for a given system, we can animate the integration process as follows:


```julia
using Mavi.Visualization

# Creating the system
system = ...

animate(system, step!)
```

It is possible to configure aspects of the animation by passing an instance of `AnimationCfg` to animate. The following example animates the system with the FPS set to 30, executing 15 time steps per animation frame:

```julia
using Mavi.Visualization

# Creating the system
system = ...

anim_cfg = AnimationCfg(
    fps = 30,
    num_steps_per_frame = 15,
)

animate(system, step!, anim_cfg)
```

## Extending the Information Panel
It is possible to inject custom information into the information panel. This is done by setting the `custom_items` field of `DefaultInfoUICfg`, which in turn is a field of `AnimationCfg`. `custom_items` is a function that should return the additional information to be displayed in the information panel. For more details about its signature, see the documentation in [`info_ui.jl`](src/gui/info_ui.jl).

The following example uses a system already defined in `Mavi.jl` and adds the position of the first particle to the information panel.

```julia
using Mavi
using Mavi.Configs
using Mavi.Visualization
using Printf

system = System(
    state = State{Float64}(
        pos=[[1 2 3]; [1 2 3]],
        vel=[[1 1 0]; [-1 0 2]],
    ),
    space_cfg=RectangleCfg(length=4, height=4),
    dynamic_cfg=LenJonesCfg(sigma=1, epsilon=0.1),
    int_cfg=IntCfg(dt=0.01),
)

function get_pos(system, _)
    pos = system.state.pos[:, 1]
    pos_formatted = @sprintf("(%.3f, %.3f)", pos[1], pos[2])
    return [("pos_1", pos_formatted)]
end

anim_cfg = AnimationCfg(
    info_cfg=DefaultInfoUICfg(
        custom_items=get_pos
    )
)

animate(system, Mavi.Integration.step!, anim_cfg)
```
