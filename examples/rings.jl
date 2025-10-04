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

using Mavi.Visualization.RingsGraphs

function create_system(;num_particles, num_cols, num_rows)
    interaction_cfg = Configs.HarmTruncCfg(
        k_rep=20,
        k_atr=4,
        dist_eq=1,
        dist_max=1 + 0.2,
    )

    dynamic_cfg = Configs.RingsCfg(
        p0=3.5,
        relax_time=1.0,
        vo=1.0,
        mobility=1.0,
        rot_diff=0.05,
        k_area=1.0,
        k_spring=20.0,
        l_spring=1.0,
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

    bbox = space_cfg.geometry_cfg
    max_size = interaction_cfg.dist_max * 1.1
    num_x_chunks = floor(Int, bbox.length / max_size)
    num_y_chunks = floor(Int, bbox.height / max_size)

    system = RingsSystem(
        state=state,
        space_cfg=space_cfg,
        dynamic_cfg=dynamic_cfg,
        int_cfg=Configs.RingsIntCfg(
            dt=0.01,
            p_chunks_cfg=Configs.ChunksCfg(num_x_chunks, num_y_chunks),
            device=Configs.Threaded(),
        ),
    )

    return system
end

function main(test=false)
    system = create_system(
        num_particles=10,
        num_cols=13,
        num_rows=13,
    )

    anim_cfg = AnimationCfg(
        num_steps_per_frame=15,
        graph_cfg=CircleGraphCfg(),
    )
    
    if !test
        animate(system, anim_cfg)
    else
        Rings.Mavi.run(system, tf=1)
    end
end

end

if !((@isdefined TEST_EX) && TEST_EX)
    Example.main()
end