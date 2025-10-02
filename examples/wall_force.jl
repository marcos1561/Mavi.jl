"""
Example: Wall with forces
"""
module Example

using Mavi
using Mavi.Configs
using Mavi.InitStates
using Mavi.Visualization

using StaticArrays, StatsBase

function main(test=false)
    # Init state configs
    num_p_x = 10
    num_p_y = 10
    offset = 0.4

    num_p = num_p_x * num_p_y
    dynamic_cfg = HarmTruncCfg(10, 1, 1)
    radius = particle_radius(dynamic_cfg)

    pos, geometry_cfg = rectangular_grid(num_p_x, num_p_y, offset, radius)
    state = Mavi.SecondLawState(
        pos=pos,
        vel=random_vel(num_p, 1/5),
    )

    geometry_cfg = geometry_cfg + RectangleCfg(
        length=geometry_cfg.length,
        height=geometry_cfg.height,
        bottom_left=geometry_cfg.bottom_left + SVector(geometry_cfg.length, 0)
    )

    circle_cfg = CircleCfg(
        radius=3*radius,
        center=geometry_cfg.size/2 + SVector(geometry_cfg.length/4, 0),
    )

    l, h = geometry_cfg.size
    lines_cfg = LinesCfg([
        [SVector(3/4 * l, 1/4*h), SVector(3/4 * l, 3/4*h)],
        [SVector(1/2 * l, 1/2*h), SVector(3/4 * l, 1/2*h)],
    ])
    
    chunks_cfg = ChunksCfg(
        num_cols=trunc(Int, 2*num_p_x*0.9), 
        num_rows=trunc(Int, 2*num_p_y*0.9),
    )
    
    int_cfg = IntCfg(
        dt=0.001,
        chunks_cfg=chunks_cfg,
    )

    system = System(
        state=state, 
        space_cfg=SpaceCfg((
            (RigidWalls(), geometry_cfg),
            (PotentialWalls(HarmTruncCfg(20.0, radius, radius)), circle_cfg),
            (PotentialWalls(HarmTruncCfg(20.0, radius, radius)), lines_cfg),
        )),
        dynamic_cfg=dynamic_cfg,
        int_cfg=int_cfg,
    )

    anim_cfg = Visualization.AnimationCfg(
        num_steps_per_frame=200,
        exec_times_size=100,
        graph_cfg=MainGraphCfg((
            CircleGraphCfg(colors_map=:magma),
            # NumsGraphCfg(),
        ))
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