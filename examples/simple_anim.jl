"""
Creates a system with two particles and animate it using 
a step function with constant velocity and rigid walls. 
"""

using Mavi: Mavi, Visualization
using Mavi.Configs

state = Mavi.State{Float64}(
    x=[1, 2],
    y=[1, 2],
    vx=[0.3, 2],
    vy=[1, 0],
)

system = Mavi.System(
    state=state, 
    space_cfg=RectangleCfg(
        length=5,
        height=5,
    ), 
    dynamic_cfg=HarmTruncCfg(1, 1, 1),
    int_cfg=IntegrationCfg(
        dt=0.1
    ),
)

function step!(system::Mavi.System)
    state = system.state
    dt = system.int_cfg.dt

    # Constant velocity
    state.x .+= state.vx * dt 
    state.y .+= state.vy * dt 

    # Rigid walls collisions
    space_cfg = system.space_cfg
    r = system.dynamic_cfg.ro/2
    for i in 1:system.num_p
        if ((state.x[i]+r) > space_cfg.length) || ((state.x[i]-r) < 0)
            state.vx[i] *= -1.
        end
        if ((state.y[i]+r) > space_cfg.height) || ((state.y[i]-r) < 0)
            state.vy[i] *= -1.
        end
    end
end

anim_cfg = Visualization.AnimationCfg(
    fps = 60,    
    circle_radius = particle_radius(system.dynamic_cfg),
    num_steps_per_frame = 1,
)

Visualization.animate(system, step!, anim_cfg)