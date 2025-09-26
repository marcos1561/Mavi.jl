module TestRings

using Random

using Mavi.Configs
using Mavi.Rings.Integration: step!

include("rings_utils.jl")
using .RingsUtils

function check_chunks(system_name)
    rng1 = MersenneTwister(24042001)
    rng2 = MersenneTwister(24042001)

    s1 = create_systems_list[system_name](rng1, Sequencial(), true)
    s2 = create_systems_list[system_name](rng2, Sequencial(), false)

    while s1.time_info.time < 10
        step!(s1)
        step!(s2)
    end

    state_square_distance(s1, s2)
end

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

function run_tests()
    # check_chunks(:normal, Threaded())
    check_threaded(:types)
end

end # TestRings

# import .TestRings
# TestRings.run_tests()