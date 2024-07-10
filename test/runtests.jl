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
            space_cfg=RectangleCfg(1, 1), 
            dynamic_cfg=HarmTruncCfg(1, 1, 1),
            int_cfg=IntegrationCfg(1)
        )
        return
    end

    function test_space_check(;state, space_cfg)
        system = Mavi.System(
            state=state, 
            space_cfg=space_cfg, 
            dynamic_cfg=DynamicCfg(1, 1, 1),
            int_cfg=IntegrationCfg(1)
        )
    end

    @test create_system() === nothing

    # @test_throws test_space_check(
    #     state = Mavi.State{Float64}(
    #         x=[1, 2, 2.5],
    #         y=[1, 2, 6],
    #         vx=[1, 2, 0],
    #         vy=[1, 2, 0],
    #     ),
    #     space_cfg = RectangleCfg(3, 5)
    # )
end
