"""
Example creating many particles to demonstrate the difference
in performance between using chunks or not.
"""

using Mavi: Mavi, Visualization
using Mavi.Configs
using Mavi.InitStates

function main()
    # Set to false to see the performance without chunks.
    use_chunks = true

    # Init state configs
    num_p_x = 25
    num_p_y = 25
    offset = 0.4

    num_p = num_p_x * num_p_y
    dynamic_cfg = HarmTruncCfg(10, 1, 1)
    radius = particle_radius(dynamic_cfg)

    pos, geometry_cfg = InitStates.rectangular_grid(num_p_x, num_p_y, offset, radius)
    state = Mavi.SecondLawState{Float64}(
        pos=pos,
        vel=InitStates.random_vel(num_p, 1/5),
    )

    if use_chunks
        int_cfg = ChunksIntCfg(
            dt=0.001,
            chunks_cfg=ChunksCfg(
                num_cols=trunc(Int, num_p_x*0.9), 
                num_rows=trunc(Int, num_p_y*0.9),
            ),
        ) 
    else
        int_cfg = int_cfg=IntCfg(dt=0.001)
    end


    system = Mavi.System(
        state=state, 
        space_cfg=SpaceCfg(
            wall_type=RigidWalls(),
            geometry_cfg=geometry_cfg,
        ),
        dynamic_cfg=dynamic_cfg,
        int_cfg=int_cfg,
    )

    anim_cfg = Visualization.AnimationCfg(
        fps=60,    
        num_steps_per_frame=300,
        exec_times_size=100,
    )

    Visualization.animate(system, Mavi.Integration.step!, anim_cfg)
end
main()