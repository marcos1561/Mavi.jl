module SerderTests

using StaticArrays, Random

using Mavi.Rings
using Mavi.Rings.Sources
using Mavi.Rings.States
using Mavi.Rings.Configs
using Mavi.Rings.InitStates

using Mavi.Visualization

include("rings_utils.jl")
using .RingsUtils

function serder_rings(name)
    root = joinpath(@__DIR__, "../tmp/serder_rings_test_$name")

    isdir(root) && rm(root; recursive=true)

    system = create_systems_list[name](MersenneTwister(24042001))

    for _ in 1:100
        Rings.Integration.step!(system)
    end
    
    save_system(system, root)
    load_system(root)
    
    for _ in 1:100
        Rings.Integration.step!(system)
    end
    
    rm(root; recursive=true)
    
    return
end

"""
Test system serialization and deserialization. Check if, after deserialization,
the system keeps evolving as if nothing has happen. 
"""
function ring_reproducibility_test(name)
    function run_sim(; serder, name)
        system = create_systems_list[name](MersenneTwister(24042001))
    
        while system.time_info.time < 10
            Rings.Integration.step!(system)
        end
        
        root = joinpath(@__DIR__, "../tmp/serder_rings_reproducibility_test_$name")

        isdir(root) && rm(root; recursive=true)
        
        if serder
            save_system(system, root)
            system = load_system(root)
        end

        while system.time_info.time < 20
            Rings.Integration.step!(system)
        end
        
        isdir(root) && rm(root; recursive=true)
        
        return system
    end

    s1 = run_sim(name=name, serder=false)
    s2 = run_sim(name=name, serder=true)
    return state_square_distance(s1, s2)
end

end # SerderTests

# import .SerderTests
# SerderTests.serder_ring()
# SerderTests.ring_reproducibility_test()
