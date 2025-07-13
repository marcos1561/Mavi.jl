"""
Simple exemple using Run-and-Tumble model
"""

using Mavi
using Mavi.States
using Mavi.InitStates
using Mavi.Configs
using Mavi.Visualization

function main()
    num_particles = 50

    num_p_x = trunc(Int, sqrt(num_particles))
    num_p_y = trunc(Int, sqrt(num_particles))

    dynamic_cfg = RunTumbleCfg(
        vo=1,
        sigma=1,
        epsilon=1,
        tumble_rate=1,
    )

    pos, geometry_cfg = rectangular_grid(
        num_p_x,
        num_p_y,
        1, # offset
        particle_radius(dynamic_cfg), # radius
    )

    space_cfg = SpaceCfg(
        wall_type=PeriodicWalls(),
        geometry_cfg=geometry_cfg,
    )

    system = System(
        state=SelfPropelledState(
            pos=pos,
            pol_angle=rand(Float64, size(pos)[2])*2*π,
        ),
        space_cfg=space_cfg,
        dynamic_cfg=dynamic_cfg,
        int_cfg=ChunksIntCfg(
            dt=0.001,
            chunks_cfg=ChunksCfg(
                num_cols=num_p_x-1, 
                num_rows=num_p_y-1, 
            ),
        ),
    )

    animate(
        system, 
        Mavi.Integration.rtp_step!, 
        AnimationCfg(num_steps_per_frame=100),
    )
end
main()