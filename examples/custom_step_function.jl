"""
Example: Custom step function

Creates a system with two particles and animate it using 
a step function with constant velocity and rigid walls. 
"""
module Example

using Mavi.StaticArrays

using Mavi
using Mavi.States
using Mavi.Visualization
using Mavi.Configs


function main(test=false)
    state = SecondLawState(
        pos=[1 2; 1 2],
        vel=[0.3 2; 1 0]
    )

    system = System(
        state=state, 
        space_cfg=SpaceCfg(
            wall_type=RigidWalls(),
            geometry_cfg=RectangleCfg(
                length=5,
                height=5,
            ), 
        ),
        dynamic_cfg=HarmTruncCfg(1, 1, 1),
        int_cfg=IntCfg(
            dt=0.1
        ),
    )

    function my_step!(system::System)
        state = system.state
        dt = system.int_cfg.dt

        # Constant velocity
        state.pos .+= state.vel * dt

        # Rigid walls collisions
        geometry_cfg = system.space_cfg.geometry_cfg
        r = system.dynamic_cfg.ro/2
        for i in 1:get_num_total_particles(system)
            if ((state.pos[i].x+r) > geometry_cfg.length) || ((state.pos[i].x-r) < 0)
                state.vel[i] = state.vel[i] .* SVector(-1, 1)
            end
            if ((state.pos[i].y+r) > geometry_cfg.height) || ((state.pos[i].y-r) < 0)
                state.vel[i] = state.vel[i] .* SVector(1, -1)
            end
        end
    end

    anim_cfg = AnimationCfg(
        num_steps_per_frame=1,
        fps=60,
    )
    if !test
        animate(system, anim_cfg, step_func=my_step!)
    else
        Mavi.run_system(system, num_steps=100, step_func=my_step!)
    end
end

end

if !((@isdefined TEST_EX) && TEST_EX)
    Example.main()
end