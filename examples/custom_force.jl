"""
Example: Custom Force

This example demonstrates how to define and use a custom pairwise force. 
A new force type, `RadialForce`, is implemented, which applies a constant 
repulsive force between particles within a specified minimum distance. 
The simulation initializes a grid of particles and visualizes their 
dynamics under this custom interaction.
"""
module Example

using Mavi
using Mavi.Configs
using Mavi.Integration
using Mavi.Visualization

struct RadialForce <: DynamicCfg
    force::Float64
    min_dist::Float64
end

Configs.particle_radius(dynamic_cfg::DynamicCfg) = dynamic_cfg.min_dist/2

function Integration.calc_interaction(i, j, dynamic_cfg::RadialForce, system)
    pos = system.state.pos
    dr = calc_diff(pos[i], pos[j], system.space_cfg)
    dist = sqrt(sum(dr.^2))

    force = dynamic_cfg.force
    min_dist = dynamic_cfg.min_dist
    
    if dist > min_dist
        return zero(dr)
    end

    return force * dr / dist
end

function main()
    # Init state configs
    num_p_x = 10
    num_p_y = 10
    offset = 0.4

    num_p = num_p_x * num_p_y
    dynamic_cfg = RadialForce(1, 1)
    radius = particle_radius(dynamic_cfg)

    pos, geometry_cfg = Mavi.InitStates.rectangular_grid(num_p_x, num_p_y, offset, radius)
    state = Mavi.SecondLawState(
        pos=pos,
        vel=Mavi.InitStates.random_vel(num_p, 1/5),
    )

    int_cfg = IntCfg(
        dt=0.001,
        chunks_cfg=ChunksCfg(
            num_cols=trunc(Int, num_p_x*0.9), 
            num_rows=trunc(Int, num_p_y*0.9),
        ),
    )

    system = System(
        state=state, 
        space_cfg=SpaceCfg(
            wall_type=RigidWalls(),
            geometry_cfg=geometry_cfg,
        ),
        dynamic_cfg=dynamic_cfg,
        int_cfg=int_cfg,
    )

    anim_cfg = AnimationCfg(
        num_steps_per_frame=300,
        graph_cfg=CircleGraphCfg(colors_map=:viridis),
    )

    animate(system, anim_cfg)
end

end

import .Example
Example.main()