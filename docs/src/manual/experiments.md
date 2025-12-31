# Running Experiments
To collect data with simulations, we need to run an experiment. To create an experiment, one needs a system and a `Collector` (object that collects data while the simulation is running). Let's demonstrate how to run an experiment in `Mavi`:

Given we have a function `create_system` (which creates a system), the first thing to do is to create the experiment configurations

```julia
using Mavi.Experiments

system = create_system()

experiment_cfg = ExperimentCfg(
    tf=10, # Final simulation time
    root="experiment_data", # Where the experiment data will be saved
)
```

next, we need to create the collector configurations, it's possible to use existing ones or create a custom collector, let's use the `DelayedCol`, which always have the system state some time in the past (this collector is useful to have the state system just before something bad happens, allowing the coder to investigate what is happening)

```julia
using Mavi.Experiments

system = create_system()

experiment_cfg = ExperimentCfg(
    tf=10, # Final simulation time
    root="experiment_data", # Where the experiment data will be saved
)

collector_cfg = DelayedCfg(
    delay_time=4, # How far in the past the state is
)
```

now we can create an experiment and run it

```julia
using Mavi.Experiments

system = create_system()

experiment_cfg = ExperimentCfg(
    tf=10, # Final simulation time
    root="experiment_data", # Where the experiment data will be saved
)

collector_cfg = DelayedCfg(
    delay_time=4, # How far in the past the state is
)

experiment = Experiment(experiment_cfg, collector_cfg, system)

run_experiment(experiment)
```

all data collected by the `Collector` and other things will be inside in the `root` path, specified when creating `ExperimentCfg`. Here's a list (not complete) of what is saved in an experiment:

- data: data collected by the `Collector`.
- final_state: system state just after the experiment finished.
- configs.json: system configurations used.
- experiment_configs.json: experiment configurations used.

## Using checkpoints
If your experiment takes a lot of time, it's a good ideia to run the experiment with checkpoints, in that way, an experiment can be resumed from a checkpoint if something bad happens. To do that, simply pass a `CheckpointCfg` to `ExperimentCfg`

```julia
experiment_cfg = ExperimentCfg(
    tf=10,
    root="experiment_data",
    checkpoint_cfg=CheckpointCfg(delta_time=2), # Creates a checkpoint every two seconds (in simulation time)
)
```

and here is how to resume an experiment from a checkpoint

```julia
using Mavi.Experiment

experiment = load_experiment("path-to-experiment")
run_experiment(experiment)
```
