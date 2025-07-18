using Mavi
using Mavi.Configs
using Test

@testset "Mavi.jl" begin
    function create_system()
        state = Mavi.SecondLawState(
            pos=[[1, 1];; [2, 1];; [3, 2]],
            vel=[[0, 0];; [1, 1];; [-1, 2]],
        )

        system = System(
            state=state, 
            space_cfg=SpaceCfg(
                geometry_cfg=RectangleCfg(length=5, height=5),
                wall_type=RigidWalls(),
            ),
            dynamic_cfg=HarmTruncCfg(1, 1, 1),
            int_cfg=IntCfg(dt=1),
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
