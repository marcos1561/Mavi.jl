using Mavi
using Mavi.Rings
using Mavi.Configs
using Test

include("test_examples.jl")

@testset "Mavi.jl" begin
    @test Examples.test()
end

include("test_rings.jl")

@testset "Rings" begin
    @test TestRings.check_chunks(:normal) < 1e-4
    @test TestRings.check_chunks(:source) < 1e-4
    @test TestRings.check_chunks(:types) < 1e-4
    
    @test TestRings.check_threaded(:normal) < 1e-4
    @test TestRings.check_threaded(:source) < 1e-4
    @test TestRings.check_threaded(:types) < 1e-4
    
    @test TestRings.check_particle_neighbors_run(:normal) === nothing
    @test TestRings.check_particle_neighbors_run(:source) === nothing
    @test TestRings.check_particle_neighbors_run(:types) === nothing
    
    @test TestRings.check_particle_neighbors_all(
        type=:rings,
        expected=[
            (pid=158, list=[142, 143]),
            (pid=222, list=[13]),
            (pid=168, list=Int64[]),
            (pid=202, list=[20, 6, 7]),
        ],
    ) == true
    @test TestRings.check_particle_neighbors_all(
        type=:all,
        expected=[
            (pid=158, list=[157, 142, 143, 159]),
            (pid=222, list=[223, 221, 13]),
            (pid=168, list=[167, 169]),
            (pid=202, list=[20, 6, 203, 201, 7]),
        ],
    ) == true
end