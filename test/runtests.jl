using Mavi
using Mavi.Configs
using Test

@testset "Mavi.jl" begin
    function create_system()
        state = Mavi.State{Float64}(
            x=[1,2],
            y=[1,2],
            vx=[1,2],
            vy=[1,2],
        )

        system = Mavi.System(
            state=state, 
            space_cfg=SpaceCfg(1, 1), 
            dynamic_cfg=DynamicCfg(1, 1, 1),
            int_cfg=IntegrationCfg(1)
        )
        return
    end

    @test create_system() === nothing
end
