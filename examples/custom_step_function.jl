"""
Example: Custom step function

Creates a system with two particles and animate it using 
a step function with constant velocity and rigid walls. 
"""
module Example

using Mavi
using Mavi.Visualization
using Mavi.Configs

function main()
    state = Mavi.States.SecondLawState{Float64}(
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

    function step!(system::System)
        state = system.state
        dt = system.int_cfg.dt

        # Constant velocity
        state.pos .+= state.vel * dt

        # Rigid walls collisions
        geometry_cfg = system.space_cfg.geometry_cfg
        r = system.dynamic_cfg.ro/2
        for i in 1:system.num_p
            if ((state.pos[1, i]+r) > geometry_cfg.length) || ((state.pos[1, i]-r) < 0)
                state.vel[1, i] *= -1.
            end
            if ((state.pos[2, i]+r) > geometry_cfg.height) || ((state.pos[2, i]-r) < 0)
                state.vel[2, i] *= -1.
            end
        end
    end

    anim_cfg = AnimationCfg(
        fps = 60,    
        num_steps_per_frame = 1,
    )

    animate(system, step!, anim_cfg)
end

end

import .Example
Example.main()