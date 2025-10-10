"""
Example: Running an Experiment

This example shows how to run an experiment to collect data from simulations.
"""
module Example

using Mavi
using Mavi.States
using Mavi.InitStates
using Mavi.Configs
using Mavi.Visualization
using Mavi.Experiments
using Mavi.Utils.Progress

function create_system()
    num_particles = 100

    num_p_x = trunc(Int, sqrt(num_particles))
    num_p_y = trunc(Int, sqrt(num_particles))

    dynamic_cfg = SzaboCfg(
        vo=1.,
        mobility=1.,
        relax_time=1.,
        k_rep=10.,
        k_adh=0.75,
        r_eq=1.,
        r_max=1+0.1,
        rot_diff=0.01,
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

    System(
        state=SelfPropelledState(
            pos=pos,
            pol_angle=rand(Float64, length(pos))*2*π,
        ),
        space_cfg=space_cfg,
        dynamic_cfg=dynamic_cfg,
        int_cfg=IntCfg(
            dt=0.01,
            chunks_cfg=ChunksCfg(
                num_cols=num_p_x-1, 
                num_rows=num_p_y-1,
            ),
        ),
    )
end

function main(test=false; delete_data=false, verbose=true)
    system = create_system()

    experiment_cfg = ExperimentCfg(
        tf=10,
        root="example_experiment_data",
        checkpoint_cfg=CheckpointCfg(delta_time=2),
    )

    collector_cfg = DelayedCfg(
        delay_time=4,
    )

    experiment = Experiment(experiment_cfg, collector_cfg, system)

    formatter = verbose ? NormalFormatter() : SilenceFormatter()
    run_experiment(experiment, prog_kwargs=(formatter=formatter,))

    if delete_data
        isdir(experiment_cfg.root) && rm(experiment_cfg.root, recursive=true)
    end

    # Discommend to animate the system
    # animate(system)
end

end #Tmp

if !((@isdefined TEST_EX) && TEST_EX)
    Example.main()
end