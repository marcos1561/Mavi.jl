module TestRings

using Random

using Mavi.Configs
using Mavi.Rings.NeighborsMod
using Mavi.Rings.Integration: step!

using Mavi.Visualization
# using Mavi.Visualization.RingsGraphs

include("rings_utils.jl")
using .RingsUtils

"Checks if integration is the same with and without chunks."
function check_chunks(system_name)
    rng1 = MersenneTwister(24042001)
    rng2 = MersenneTwister(24042001)

    s1 = create_systems_list[system_name](rng1, Sequencial(), use_chunks=true)
    s2 = create_systems_list[system_name](rng2, Sequencial(), use_chunks=false)

    while s1.time_info.time < 10
        step!(s1)
        step!(s2)
    end

    state_square_distance(s1, s2)
end

"Checks if integration is the same with and without threads."
function check_threaded(system_name)
    rng1 = MersenneTwister(24042001)
    rng2 = MersenneTwister(24042001)

    s1 = create_systems_list[system_name](rng1, Threaded())
    s2 = create_systems_list[system_name](rng2, Sequencial())

    while s1.time_info.time < 10
        step!(s1)
        step!(s2)
    end

    state_square_distance(s1, s2)
end

function check_particle_neighbors_run(system_name)
    rng = MersenneTwister(1532466)
    s = create_systems_list[system_name](rng, Threaded(),
        p_neighbors_cfg=NeighborsCfg(),
    )

    while s.time_info.time < 4
        step!(s)
    end

    return
end

function check_particle_neighbors(device, cfg_dict; use_chunks, expected, show=false)
    rng = MersenneTwister(31415)

    s = create_system_normal(
        rng, device,
        use_chunks=use_chunks,
        num_particles=10,
        num_cols=5,
        num_rows=5;
        cfg_dict...,
    )

    cfg = cfg_dict[1]

    while s.time_info.time < 10
        step!(s)
    end

    if :p_neighbors_cfg in keys(cfg_dict)
        neigh = s.info.p_neigh
    else
        neigh = s.info.r_neigh
    end

    if !show
        for pid_exp in expected
            pid = pid_exp.pid
            list_exp = sort(pid_exp.list)
            num_exp = length(list_exp)

            num_actual = get_neigh_count(neigh)[pid]

            if cfg.only_count
                list_actual = get_neigh(neigh).list
                list_exp = nothing
            else
                list_actual = sort(get_neigh_list(neigh, pid))
            end

            if (num_exp != num_actual) || (list_exp != list_actual)
                return false
            end
        end
        return true
    end

    if show
        pids = [158, 222, 168, 202]
        # pids = [8, 20, 11]    

        for pid in pids
            println("\n$(pid):")
            println("Num :", get_neigh_count(neigh)[pid])
            println("List:", get_neigh_list(neigh, pid))
        end

        anim_cfg = AnimationCfg(
            graph_cfg=MainGraphCfg((
                CircleGraphCfg(rng=MersenneTwister(12345678)),
                RingsNumsGraphCfg(),
            )),
            begin_paused=true,
        )
        animate(s, anim_cfg)
    end
end

"""
Checks particles neighbors calculation for all possible devices
and chunks method (with/without chunks). 

For every possible configuration (device, use_chunks) does the fallowing:
create a system with the given configurations, run until some
final time and checks if the neighbors are equal to `expected`.
    
# Arguments
- `type`:  
    Neighbors calculation type.
- `expected`:  
    List of `NamedTuple` with two fields:  
    - pid: Particle index.
    - list: List of expected neighbors indices after the
        simulation is run.

# Return
`true` if no erros, `false` otherwise.
"""
function check_particle_neighbors_all(; type, expected)
    cfg_only_count = (p_neighbors_cfg=NeighborsCfg(
        only_count=true,
        type=type,
        tol=1.1,
    ),)
    cfg_not_only_count = (p_neighbors_cfg=NeighborsCfg(
        only_count=false,
        type=type,
        tol=1.1,
    ),)

    devices = [Sequencial(), Threaded()]
    use_chunks = [true, false]

    for device in devices
        for chunk in use_chunks
            if device == Threaded() && !chunk
                continue
            end

            if !check_particle_neighbors(device, cfg_only_count, use_chunks=chunk, expected=expected)
                println("Error: device=$device | use_chunks=$chunk | only_count=true")
                return false
            end
            if !check_particle_neighbors(device, cfg_not_only_count, use_chunks=chunk, expected=expected)
                println("Error: device=$device | use_chunks=$chunk | only_count=false")
                return false
            end
        end
    end

    return true
end

function run_tests()
    # check_chunks(:normal, Threaded())
    # check_threaded(:types)

    # t1 = check_particle_neighbors_all(
    #     type=:rings,
    #     expected=[
    #         (pid=158, list=[142, 143]),
    #         (pid=222, list=[13]),
    #         (pid=168, list=Int64[]),
    #         (pid=202, list=[20, 6, 7]),
    #     ],
    # )
    # t2 = check_particle_neighbors_all(
    #     type=:all,
    #     expected=[
    #         (pid=158, list=[157, 142, 143, 159]),
    #         (pid=222, list=[223, 221, 13]),
    #         (pid=168, list=[167, 169]),
    #         (pid=202, list=[20, 6, 203, 201, 7]),
    #     ],
    # )
    # println(t1, ", ", t2)

    # check_particle_neighbors_run(:normal)
    # check_particle_neighbors_run(:source)
    # check_particle_neighbors_run(:types)

    # run_examples()
    Examples.test()

    # check_particle_neighbors(
    #     Threaded(),
    #     (r_neighbors_cfg=NeighborsCfg(
    #         only_count=false,
    #         tol=1.1,
    #     ),),
    #     use_chunks=true,
    #     expected=[
    #         (pid=158, num=2, list=[142, 143]),
    #         (pid=222, num=1, list=[13]),
    #         (pid=168, num=0, list=Int64[]),
    #         (pid=202, num=3, list=[20, 6, 7]),
    #     ],
    #     show=true,
    # )
end

end # TestRings

# import .TestRings
# TestRings.run_tests()