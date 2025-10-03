"""
Example: Using Chunks

This example creates many particles to demonstrate the difference
in performance between using chunks or not.
"""
module Example

using Mavi
using Mavi.Configs
using Mavi.InitStates
using Mavi.Visualization

function main(test=false)
    # Set to false to see the performance without chunks.
    use_chunks = true

    # Init state configs
    num_p_x = 25
    num_p_y = 25
    offset = 0.4

    num_p = num_p_x * num_p_y
    dynamic_cfg = HarmTruncCfg(10, 1, 1)
    radius = particle_radius(dynamic_cfg)

    pos, geometry_cfg = rectangular_grid(num_p_x, num_p_y, offset, radius)
    state = Mavi.SecondLawState(
        pos=pos,
        vel=random_vel(num_p, 1/5),
    )

    chunks_cfg = nothing
    if use_chunks    
        chunks_cfg = ChunksCfg(
            num_cols=trunc(Int, num_p_x*0.9), 
            num_rows=trunc(Int, num_p_y*0.9),
        )
    end
    
    int_cfg = IntCfg(
        dt=0.001,
        chunks_cfg=chunks_cfg,
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

    anim_cfg = Visualization.AnimationCfg(
        num_steps_per_frame=200,
        exec_times_size=100,
        graph_cfg=CircleGraphCfg(colors_map=:magma),
    )

    if !test
        animate(system, anim_cfg)
    else
        Mavi.run(system, tf=1)
    end
end

end

if !((@isdefined TEST_EX) && TEST_EX)
    Example.main()
end