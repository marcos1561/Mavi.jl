"""
Example: Rings

Simple simulation using the Rings module. Description:
- Only one ring type.
- The number of rings is constant.
- Rings are created in a regular rectangular grid with random polarizations.
"""
module Example

using Mavi.Rings
using Mavi.Visualization

function create_system(;num_particles, num_cols, num_rows)
    interaction_cfg = Configs.HarmTruncCfg(
        k_rep=13,
        k_atr=1,
        dist_eq=1,
        dist_max=1 + 0.1,
    )

    dynamic_cfg = Configs.RingsCfg(
        p0=3.5,
        relax_time=1.0,
        vo=1.0,
        mobility=1.0,
        rot_diff=0.05,
        k_area=1.0,
        k_spring=10.0,
        l_spring=1.0,
        num_particles=num_particles,
        interaction_finder=interaction_cfg,
    )

    num_rings = num_cols * num_rows
    rings_pos, geometry_cfg = InitStates.rectangular_grid(
        num_cols=num_cols,
        num_rows=num_rows,
        num_particles=num_particles,
        p_radius=Configs.particle_radius(dynamic_cfg),
        pad_x=0.1, pad_y=0.1,
    )

    space_cfg = Configs.SpaceCfg(
        geometry_cfg=geometry_cfg,
        wall_type=Configs.PeriodicWalls(),
    )
    state = RingsState(
        rings_pos=rings_pos,
        pol=InitStates.random_pol(num_rings),
    )

    system = RingsSystem(
        state=state,
        space_cfg=space_cfg,
        dynamic_cfg=dynamic_cfg,
        int_cfg=Configs.IntCfg(
            dt=0.01,
            chunks_cfg=Configs.ChunksCfg(44, 44)
        ),
    )

    print(typeof(system))

    return system
end

function main()
    system = create_system(
        num_particles=10,
        num_cols=13,
        num_rows=13,
    )
    
    colors = []
    for _ in 1:size(system.state.rings_pos, 3)
        push!(colors, rand(Visualization.RGBf))
    end
    
    anim_cfg = AnimationCfg(
        num_steps_per_frame=10,
        graph_cfg=DefaultGraphCfg(
            colors_map=colors,
        ),
    )

    animate(system, Integration.step!, anim_cfg)
end

end

import .Example
Example.main()