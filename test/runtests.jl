using Test

include("tests_examples.jl")
include("tests_experiments.jl")
@testset "Mavi Core" begin
    @test Examples.test()
    
    @testset "Experiments" begin
        @test TestExperiments.experiment_batch() == true
        @test TestExperiments.load_batch() == true
        @test TestExperiments.add_experiments_to_batch() == true
    end 
end

include("tests_rings/runtests.jl")