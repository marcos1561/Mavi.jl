module VisualizationTests

using Test
using Mavi
using Mavi.StaticArrays
using Mavi.Visualization
using Mavi.Visualization.SystemGraphs
using Mavi.Visualization: GLMakie
using GLMakie: Figure, Axis, Point2f, Vec2f, RGBf, RGBAf

const SG = Mavi.Visualization.SystemGraphs

struct TestRadii <: Mavi.Configs.DynamicCfg
    values::Vector{Float64}
end
Mavi.Systems.get_particle_radius(cfg::TestRadii, state, id) = cfg.values[id]

function create_system(; active_state=nothing)
    state = Mavi.States.SecondLawState(
        pos=[SVector(1.0, 1.0), SVector(3.0, 2.0), SVector(5.0, 3.0)],
        vel=fill(SVector(0.0, 0.0), 3),
        active_state=active_state,
    )
    Mavi.Systems.System(
        state=state,
        dynamic_cfg=TestRadii([0.25, 0.5, 0.75]),
        space_cfg=Mavi.Configs.SpaceCfg(
            geometry_cfg=Mavi.Configs.RectangleCfg(length=8.0, height=8.0),
            wall_type=Mavi.Configs.PeriodicWalls(),
        ),
        int_cfg=Mavi.Configs.IntCfg(dt=0.01),
    )
end

function check_circle_graph()
    palette = RGBf[RGBf(1, 0, 0), RGBf(0, 1, 0), RGBf(0, 0, 1)]
    expected_colors = RGBAf.(palette)
    painter_cfg = PalettePainterCfg(palette=palette)

    @testset "Circle sizes, colors and in-place motion" begin
        system = create_system()
        fig = Figure()
        graph = get_graph(Axis(fig[1, 1]), system, MainGraphCfg(
            CircleGraphCfg(painter_cfg=painter_cfg),
        ))
        comp = only(graph.components)
        @test comp.plot[1][] == Point2f.(system.state.pos)
        @test comp.plot.markersize[] == Vec2f[Vec2f(0.5), Vec2f(1.0), Vec2f(1.5)]
        @test comp.plot.color[] == expected_colors

        # Check that reusing an input array still notifies Makie's downstream data.
        changes = Ref(0)
        handle = GLMakie.on(comp.plot[1]) do _
            changes[] += 1
        end
        system.state.pos[2] = SVector(4.0, 4.0)
        system.dynamic_cfg.values[2] = 0.625
        update_graph(graph, system)
        @test changes[] > 0
        @test comp.plot[1][][2] == Point2f(4, 4)
        @test comp.plot.markersize[][2] == Vec2f(1.25)
        GLMakie.off(handle)
    end

    @testset "Dynamic particle colors" begin
        system = create_system()
        dynamic_painter = PalettePainterCfg(
            palette=palette,
            types=DynamicTypes(types=:entities, updater=(painter, system) -> begin
                painter.types[1] = system.state.pos[1].x > 2 ? 3 : 1
                return true
            end),
        )
        fig = Figure()
        graph = get_graph(Axis(fig[1, 1]), system, MainGraphCfg(
            CircleGraphCfg(painter_cfg=dynamic_painter),
        ))
        comp = only(graph.components)
        @test comp.plot.color[][1] == expected_colors[1]
        system.state.pos[1] = SVector(4.0, 1.0)
        update_graph(graph, system)
        @test comp.plot.color[][1] == expected_colors[3]
    end

    @testset "Particles appear, disappear and reuse slots" begin
        system = create_system(active_state=Mavi.States.ActiveState(fill(false, 3)))
        fig = Figure()
        graph = get_graph(Axis(fig[1, 1]), system, MainGraphCfg(
            CircleGraphCfg(painter_cfg=painter_cfg),
        ))
        comp = only(graph.components)
        for ids in (Int[], [2, 3], [1], Int[], [1, 2, 3], [3])
            system.state.part_ids.mask .= false
            system.state.part_ids.mask[ids] .= true
            Mavi.States.update_ids!(system.state)
            update_graph(graph, system)
            @test comp.plot[1][] == Point2f.(system.state.pos[ids])
            @test comp.plot.markersize[] == Vec2f.(2 .* system.dynamic_cfg.values[ids])
            @test comp.plot.color[] == expected_colors[ids]
        end
    end

    @testset "Fixed radius and outline options" begin
        system = create_system(active_state=Mavi.States.ActiveState([false, true, true]))
        fig = Figure()
        graph = get_graph(Axis(fig[1, 1]), system, MainGraphCfg((
            CircleGraphCfg(circle_radius=0.4, stroke_cfg=nothing, painter_cfg=painter_cfg),
            CircleGraphCfg(stroke_cfg=SG.StrokeCfg(color=:red, width=2.0), painter_cfg=painter_cfg),
        )))
        fixed, outlined = graph.components
        @test fixed.plot.markersize[] == [Vec2f(0.8), Vec2f(0.8)]
        @test fixed.plot.strokewidth[] == 0
        @test outlined.plot.strokewidth[] == 2
        @test GLMakie.to_color(outlined.plot.strokecolor[]) == GLMakie.to_color(:red)
    end
end

@testset "CircleGraph rendering" begin
    check_circle_graph()
end

end
