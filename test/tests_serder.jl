using Test
import Mavi
using Mavi.MaviSerder: JSON3, register_type_module!

module SerderTestTypes
    import Mavi

    struct Config{T}
        value::T
    end

    callback() = nothing
    struct CallbackConfig{F} end
    struct SystemType <: Mavi.Systems.SystemType end
end

@testset "External modules in saved types" begin
    serder = Mavi.MaviSerder
    info = (type="SerderTestTypes.Config{Float64}",)

    @test_throws UndefVarError serder.get_saved_type(info)
    @test register_type_module!(SerderTestTypes) === SerderTestTypes
    @test register_type_module!(SerderTestTypes) === SerderTestTypes
    @test_throws ArgumentError register_type_module!(Module(:SerderTestTypes))
    @test_throws ArgumentError register_type_module!(Module(:Mavi))
    @test !isdefined(serder, :SerderTestTypes)

    @test serder.get_saved_type(info) === SerderTestTypes.Config{Float64}
    @test serder.get_saved_type((type="Tuple{SerderTestTypes.Config{Float64}, Nothing}",)) ===
        Tuple{SerderTestTypes.Config{Float64}, Nothing}
    @test serder.get_saved_type((type="SerderTestTypes.CallbackConfig{typeof(SerderTestTypes.callback)}",)) ===
        SerderTestTypes.CallbackConfig{typeof(SerderTestTypes.callback)}
    @test serder.get_saved_type((type="Mavi.Configs.RectangleCfg{2, Float64}",)) ===
        Mavi.Configs.RectangleCfg{2, Float64}
    @test_throws UndefVarError serder.get_saved_type((type="UnregisteredSerderModule.Config",))

    mktempdir() do root
        cfg = SerderTestTypes.Config(2.5)
        serder.save_component_json(cfg, root, "custom")
        restored = serder.load_component_json(joinpath(root, "custom.json"))
        @test restored isa typeof(cfg)
        @test restored.value == cfg.value

        serder.save_component_serial(cfg, root, "custom_serial")
        restored_serial = serder.load_component_serial(joinpath(root, "custom_serial"))
        @test restored_serial isa typeof(cfg)
        @test restored_serial.value == cfg.value

        configs_path = joinpath(root, "configs.json")
        open(configs_path, "w") do io
            JSON3.write(io, (
                sys_type="SerderTestTypes.SystemType()",
                custom=serder.get_obj_save_data_json(cfg),
            ))
        end
        configs = serder.load_system_configs(configs_path)
        @test configs[:sys_type] isa SerderTestTypes.SystemType
        @test configs[:custom].value == cfg.value
    end
end
