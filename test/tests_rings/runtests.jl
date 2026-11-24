include("tests_general.jl")
include("tests_serder.jl")
include("tests_sources.jl")
include("tests_neighbors.jl")
include("tests_states.jl")

@testset "Rings" begin
    @testset "Valid particle pairs" begin
        RingStateTests.check_valid_pairs()
    end

    @testset "Repeatable neighbors" begin
        NeighborTests.check_repeatable_neighbors()
        NeighborTests.check_ring_neighbors()
    end

    @testset "Source boundary chunks" begin
        SourceTests.check_boundary_chunks()
    end

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

    @testset "Serder" begin
        SerderTests.check_source_components()

        @test SerderTests.serder_rings(:normal) === nothing
        @test SerderTests.serder_rings(:source) === nothing
        @test SerderTests.serder_rings(:types) === nothing
        
        @test SerderTests.ring_reproducibility_test(:normal) < 1e-4
        @test SerderTests.ring_reproducibility_test(:source) < 1e-4
        @test SerderTests.ring_reproducibility_test(:types) < 1e-4
    end
end
