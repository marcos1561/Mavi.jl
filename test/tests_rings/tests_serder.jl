module SerderTests

import Mavi
using Test
using StaticArrays, Random
using Mavi.MaviSerder: JSON3

using Mavi.Rings
using Mavi.Rings.Sources
using Mavi.Rings.States
using Mavi.Rings.Configs
using Mavi.Rings.InitStates

using Mavi.Visualization

include("rings_utils.jl")
using .RingsUtils

struct PathComponent
    file::String
end

function Mavi.MaviSerder.load_component(::Type{PathComponent}, data::JSON3.Object, path)
    PathComponent(joinpath(dirname(path), data.file))
end

function check_source_components()
    system = create_system_source(MersenneTwister(24042001))
    source, sink = (s.cfg for s in system.info.sources)
    path = joinpath(tempdir(), "mavi_components", "configs.json")

    for configs in ((), (source,), (source, sink), (source, (sink, nothing)),
        (PathComponent("geometry.bin"),))
        encoded = Mavi.MaviSerder.get_obj_save_data_json(configs)
        data = JSON3.read(JSON3.write(encoded))
        restored = Mavi.MaviSerder.load_component(Nothing, data, path)
        @test length(restored) == length(configs)
        @test typeof(restored) == typeof(configs)
        if configs == (source, sink)
            @test restored[1].spawn_pos == source.spawn_pos
            @test restored[2].geometry_cfg == sink.geometry_cfg
        elseif configs isa Tuple{PathComponent}
            @test restored[1].file == joinpath(dirname(path), "geometry.bin")
        end
    end

    for _ in 1:3
        Rings.Integration.step!(system)
    end
    mktempdir() do root
        save_system(system, root)
        restored = load_system(root)
        @test map(s -> typeof(s.cfg), restored.info.sources) ==
            map(s -> typeof(s.cfg), system.info.sources)
        pids = collect(Mavi.States.get_particles_ids(system.state))
        @test restored.state.pos[pids] == system.state.pos[pids]
        @test restored.state.rings_ids.mask == system.state.rings_ids.mask
        @test restored.time_info.time == system.time_info.time
        for _ in 1:5
            Rings.Integration.step!(system)
            Rings.Integration.step!(restored)
        end
        pids = collect(Mavi.States.get_particles_ids(system.state))
        @test collect(Mavi.States.get_particles_ids(restored.state)) == pids
        @test restored.state.pos[pids] ≈ system.state.pos[pids]
        @test restored.state.pol ≈ system.state.pol
    end
end

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
