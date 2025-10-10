module TestExperiments

using Setfield

using Mavi
using Mavi.Configs
using Mavi.InitStates
using Mavi.Experiments
using Mavi.Visualization

function get_init_system(num_p)
    num_p_x = num_p_y = round(Int, sqrt(num_p))
    num_p = num_p_x * num_p_y

    dynamic_cfg = HarmTruncCfg(10.0, 1.0, 1.0)
    radius = particle_radius(dynamic_cfg)
    
    offset = 4 * radius

    pos, geometry_cfg = rectangular_grid(num_p_x, num_p_y, offset, radius)
    state = Mavi.SecondLawState(
        pos=pos,
        vel=random_vel(num_p, 1/5),
    )

    chunks_cfg = ChunksCfg(
        num_cols=trunc(Int, num_p_x*0.9), 
        num_rows=trunc(Int, num_p_y*0.9),
    )
    
    int_cfg = IntCfg(
        dt=0.001,
        chunks_cfg=chunks_cfg,
        device=Sequencial(),
    )

    System(
        state=state, 
        space_cfg=SpaceCfg(
            wall_type=RigidWalls(),
            geometry_cfg=geometry_cfg,
        ),
        dynamic_cfg=dynamic_cfg,
        int_cfg=int_cfg,
    )
end

function get_system(init_system, value, idx)
    system = deepcopy(init_system)
    system = @set system.dynamic_cfg.ko = Float64(value.ko)
    system = @set system.int_cfg.dt = Float64(value.dt)
    return system
end

function get_system_load_test(init_system, value, idx)
    if idx == 2
        error("Intentional error!")
    end
    get_system(init_system, value, idx)
end

"""
Run an experiment batch and checks if the expected files exists after
the simulation is done.
"""
function experiment_batch(; remove_data=true, verbose=false)
    exp_cfg = ExperimentCfg(
        tf=2000.0,
        root=joinpath(@__DIR__, "tmp/experiment_batch_run_test"),
    )

    remove_data && isdir(exp_cfg.root) && rm(exp_cfg.root, recursive=true)

    col_cfg = DelayedCfg(8)

    init_system = get_init_system(30)

    values = [(ko=ko_i, dt=dt_i) for (ko_i, dt_i) in zip(LinRange(7, 13, 3), LinRange(0.001, 0.01, 3))]

    exp_batch = ExperimentBatch(
        exp_cfg=exp_cfg,
        col_cfg=col_cfg,
        init_system=init_system,
        values=values,
    )

    run_experiment_batch(exp_batch, get_system, verbose=verbose)

    exp_datas_list = sort!(collect(readdir(joinpath(exp_cfg.root, "data"))))
    
    have_experiment_data = exp_datas_list == ["1", "2", "3"]
    have_init_system = isdir(joinpath(exp_cfg.root, "init_system"))
    have_exp_configs = isfile(joinpath(exp_cfg.root, "experiment_configs.json"))
    have_values = isfile(joinpath(exp_cfg.root, "values.bin"))

    remove_data && isdir(exp_cfg.root) && rm(exp_cfg.root, recursive=true)

    # animate(init_system)
    return have_experiment_data && have_init_system && have_exp_configs && have_values
end

"""
Run a batch with an intentional error in the middle of the batch, then reload 
the experiment and finish the batch.
"""
function load_batch(; verbose=false)
    exp_cfg = ExperimentCfg(
        tf=1,
        root=joinpath(@__DIR__, "tmp/experiment_batch_load_test"),
    )

    isdir(exp_cfg.root) && rm(exp_cfg.root, recursive=true)

    col_cfg = DelayedCfg(0.5)

    init_system = get_init_system(10)

    values = [(ko=ko_i, dt=dt_i) for (ko_i, dt_i) in zip(LinRange(7, 13, 3), LinRange(0.001, 0.01, 3))]

    exp_batch = ExperimentBatch(
        exp_cfg=exp_cfg,
        col_cfg=col_cfg,
        init_system=init_system,
        values=values,
    )

    run_experiment_batch(exp_batch, get_system_load_test, verbose=verbose)

    exp_batch = load_experiment_batch(exp_cfg.root)
    run_experiment_batch(exp_batch, get_system, verbose=verbose)

    exp_datas_list = sort!(collect(readdir(joinpath(exp_cfg.root, "data"))))

    isdir(exp_cfg.root) && rm(exp_cfg.root, recursive=true)

    return exp_datas_list == ["1", "2", "3"]
end

"""
Run a batch, then load it, add two more experiments and run it again.
"""
function add_experiments_to_batch(verbose=false)
     exp_cfg = ExperimentCfg(
        tf=1,
        root=joinpath(@__DIR__, "tmp/experiment_batch_add_test"),
    )

    isdir(exp_cfg.root) && rm(exp_cfg.root, recursive=true)

    col_cfg = DelayedCfg(0.5)

    init_system = get_init_system(10)

    values = [(ko=ko_i, dt=dt_i) for (ko_i, dt_i) in zip(LinRange(7, 13, 3), LinRange(0.001, 0.01, 3))]

    exp_batch = ExperimentBatch(
        exp_cfg=exp_cfg,
        col_cfg=col_cfg,
        init_system=init_system,
        values=values,
    )

    run_experiment_batch(exp_batch, get_system, verbose=verbose)

    exp_batch = load_experiment_batch(exp_cfg.root)
    
    new_values = [(ko=10.0, dt=0.0012), (ko=10.69, dt=0.0021)]
    exp_batch = add_experiments(exp_batch, new_values)

    run_experiment_batch(exp_batch, get_system, verbose=verbose)

    exp_datas_list = sort!(collect(readdir(joinpath(exp_cfg.root, "data"))))
    
    system4 = load_system(joinpath(exp_cfg.root, "data/4", "final_state"))
    system5 = load_system(joinpath(exp_cfg.root, "data/5", "final_state"))
    
    correct_configs_4 = system4.dynamic_cfg.ko == new_values[1].ko && system4.int_cfg.dt == new_values[1].dt
    correct_configs_5 = system5.dynamic_cfg.ko == new_values[2].ko && system5.int_cfg.dt == new_values[2].dt

    isdir(exp_cfg.root) && rm(exp_cfg.root, recursive=true)

    return exp_datas_list == ["1", "2", "3", "4", "5"] && correct_configs_4 && correct_configs_5
end

function view_final_state(root)
    for exp_root in readdir(joinpath(root, "data"); join=true)
        if !isdir(exp_root)
            continue
        end
        
        system = load_system(joinpath(exp_root, "final_state"))
        
        anim_cfg = AnimationCfg(
            begin_paused=true,
        )

        animate(system, anim_cfg)
    end
end

function run_tests()
    tests_func = [
        # experiment_batch,
        # load_batch,
        add_experiments_to_batch
    ]

    results = []

    for (i, func) in enumerate(tests_func)
        if !func()
            push!(results, "Error in test $i")
        end
    end

    println("\n=============\n")
    if length(results) == 0
        println("All tests passed!")
    else
        for mssg in results
            println(mssg)
        end
    end 
end

end # TestExperiments

# TestExperiments.run_tests()
# TestExperiments.experiment_batch()
# TestExperiments.load_batch()
# TestExperiments.add_experiments_to_batch()
# TestExperiments.load_batch()
# TestExperiments.test_add_experiments()
# TestExperiments.view_final_state("experiment_batch_test")