module SerderTests

using StaticArrays, Random

using Mavi.Rings
using Mavi.Rings.Sources
using Mavi.Rings.States
using Mavi.Rings.Configs
using Mavi.Rings.InitStates

using Mavi.Visualization

include("../rings_utils.jl")
using .RingsUtils

function serder_ring()
    rm("tmp_serder_tests"; force=true, recursive=true)
    funcs = (
        normal = create_system_normal,
        source = create_system_source,
        type = create_system_types,
    )

    for create_system in values(funcs)
        system = create_system(MersenneTwister(24042001))

        for _ in 1:100
            Rings.Integration.step!(system)
        end
        
        save_system(system, "tmp_serder_tests")
        load_system("tmp_serder_tests")
        
        for _ in 1:100
            Rings.Integration.step!(system)
        end
        
        rm("tmp_serder_tests"; force=true, recursive=true)
    end
    
    return
end

"""
Test system serialization and deserialization. Check if, after deserialization,
the system keeps evolving as if nothing has happen. 
"""
function ring_reproducibility_test()
    function run_sim(; serder, name)
        system = create_systems_list[name](MersenneTwister(24042001))
    
        while system.time_info.time < 10
            Rings.Integration.step!(system)
        end
        
        rm("tmp_serder_tests"; force=true, recursive=true)
        
        if serder
            save_system(system, "tmp_serder_tests")
            system = load_system("tmp_serder_tests")
        end

        while system.time_info.time < 20
            Rings.Integration.step!(system)
        end
        
        rm("tmp_serder_tests"; force=true, recursive=true)
        
        return system
    end

    for name in keys(create_systems_list)
        println("Testing $name...")
        s1 = run_sim(name=name, serder=false)
        s2 = run_sim(name=name, serder=true)
        @assert state_square_distance(s1, s2) < 1e-4
    end
end

end # SerderTests

import .SerderTests
# SerderTests.serder_ring()
SerderTests.ring_reproducibility_test()
