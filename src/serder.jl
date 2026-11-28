module MaviSerder

export save_system, save_system_configs, save_component_serial, save_component_json, get_obj_save_data_json
export load_system, load_dic_configs, load_component, load_component_serial, load_component_json, load_configs, load_system_configs, load_state
export register_type_module!

using Serialization, JSON3, Random, StructTypes
using StaticArrays

import Mavi
using Mavi.Systems
# using Mavi.Configs: PotentialCfg

# ==
# SerDer for Matrix
# ==

StructTypes.StructType(::Type{Matrix{T}}) where T = StructTypes.CustomStruct()

function StructTypes.lower(m::Matrix{T}) where T
    return (size=size(m), data=vec(m))
end

function StructTypes.construct(::Type{Matrix{T}}, x) where T
    elems = [isa(el, T) ? el : JSON3.read(JSON3.write(el), T) for el in x["data"]]
    reshape(elems, x["size"]...)
end

# ==
# Serialization
# ==

"""
Save `obj` in the folder `root/name` using `Serialization`. Two things are saved inside the folder:
- `type.bin`: `obj` type name 
- `data.bin`: `obj` data 
"""
function save_component_serial(obj, root, name) 
    root_path = mkpath(joinpath(root, name))
    data = get_obj_save_data_json(obj, root_path)
    serialize(joinpath(root_path, "type.bin"), data[:type])
    serialize(joinpath(root_path, "data.bin"), data[:data])
end

"Returns data to be saved inside `obj`."
get_obj_save_data(obj, path) = get_obj_save_data(obj)
get_obj_save_data(obj) = obj
# get_obj_save_data(obj::Dict) = Dict(k => get_obj_save_data_json(v) for (k, v) in obj)

"Returns data to be saved inside `obj` in the `json` format."
get_obj_save_data_json(obj, path=nothing) = (type=sprint(show, typeof(obj); context=:module => nothing), data=get_obj_save_data(obj, path))
get_obj_save_data_json(obj::Tuple, path=nothing) = [get_obj_save_data_json(x, path) for x in obj]

function save_component_json(obj, root, name)
    path = joinpath(root, "$name.json")
    open(path, "w") do io
        JSON3.pretty(io, get_obj_save_data_json(obj, path))
    end
end

function save_configs(configs, save_path)
    configs_data = get_obj_save_data(configs, dirname(save_path))
    mkpath(dirname(save_path))
    open(save_path, "w") do io
        JSON3.pretty(io, configs_data)
    end
end

function save_system_configs(system::System, root, name="configs")
    mkpath(root)
    path = joinpath(root, "$name.json")
    configs_data = Dict(
        "sys_type" => system.type,
        "space_cfg" => get_obj_save_data_json(system.space_cfg, path),
        "dynamic_cfg" => get_obj_save_data_json(system.dynamic_cfg, path),
        "int_cfg" => get_obj_save_data_json(system.int_cfg, path),
        "info" => get_obj_save_data_json(system.info, path),
        "debug_info" => get_obj_save_data_json(system.debug_info, path),
        "time_info" => get_obj_save_data_json(system.time_info, path),
    )

    open(path, "w") do io
        JSON3.pretty(io, configs_data)
    end
end

"Save a system in the path `root`."
function save_system(system::System, root; metadata=nothing)
    mkpath(root)
    if metadata !== nothing
        save_component_json(metadata, root, "metadata")
    end
    save_system_configs(system, root)
    save_component_serial(system.state, root, "state")
    serialize(joinpath(root, "state/rng.bin"), system.rng)
end


# ==
# Deserialization
# ==

const REGISTERED_TYPE_MODULES = Dict{Symbol, Module}()

"""
    register_type_module!(mod::Module)

Make `mod` available by its name when resolving saved component types and system
types. Register the root module of a package that defines custom configurations
or user data before loading a saved system. Packages can register themselves in
their `__init__` function:

```julia
function __init__()
    Mavi.MaviSerder.register_type_module!(@__MODULE__)
end
```

Registering the same module again is allowed. A different module with the same
name raises an `ArgumentError`. Registration does not import the module or add
global bindings to `MaviSerder`.
"""
function register_type_module!(mod::Module)
    name = nameof(mod)
    if haskey(REGISTERED_TYPE_MODULES, name) && REGISTERED_TYPE_MODULES[name] !== mod
        throw(ArgumentError("a different module is already registered as $name"))
    end
    if isdefined(@__MODULE__, name) && getfield(@__MODULE__, name) !== mod
        throw(ArgumentError("$name already refers to a different binding in MaviSerder"))
    end
    REGISTERED_TYPE_MODULES[name] = mod
    return mod
end

resolve_type_modules(x) = x
resolve_type_modules(name::Symbol) = haskey(REGISTERED_TYPE_MODULES, name) ? QuoteNode(REGISTERED_TYPE_MODULES[name]) : name
resolve_type_modules(expr::Expr) = Expr(expr.head, map(resolve_type_modules, expr.args)...)

eval_saved_expression(value) = eval(resolve_type_modules(Meta.parse(string(value))))

get_saved_type(load_info) = eval_saved_expression(load_info[:type])
get_saved_type(load_info::JSON3.Array) = Nothing

get_saved_data(load_info) = load_info[:data]
get_saved_data(load_info::JSON3.Array) = load_info

"Information to load an object saved using `Serialization`."
get_load_info_serial(root) = (type=deserialize(joinpath(root, "type.bin")), data=root)

load_component(T, data, path) = load_component(T, data)

"Load object of type `T` with data saved in `path`."
load_component(T::Type, path::String) = deserialize(joinpath(path, "data.bin"))

"Load object of type `T` with data in `data`."
load_component(T::Type, data::JSON3.Object) = JSON3.read(JSON3.write(data), T)

load_component(T::Type{Nothing}, _) = nothing

"Load array of objects with load information in `data`"
function load_component(T::Type{Nothing}, data::JSON3.Array, path)
    data_loaded = []
    for load_info in data
        push!(data_loaded, load_component(get_saved_type(load_info), get_saved_data(load_info), path))
    end
    return Tuple(data_loaded)
end

function load_component_serial(path)
   load_info = get_load_info_serial(path)
   load_component(get_saved_type(load_info), get_saved_data(load_info), path) 
end

function load_component_json(path)
   load_info = JSON3.read(path)
   load_component(get_saved_type(load_info), get_saved_data(load_info), path) 
end

"Load dictionary of object with load information in `configs`."
function load_dic_configs(configs, path)
    configs_loaded = Dict()
    for (name, load_info) in configs
        if name == :sys_type
            continue
        end
        # println("Name: $name")
        # println("Type: ", get_saved_type(load_info))
        # println("Data Type: ", typeof(get_saved_data(load_info)))
        # println("Data Value: ", get_saved_data(load_info))
        # println("Data: ", JSON3.write(get_saved_data(load_info)))
        # println("================")
        configs_loaded[name] = load_component(get_saved_type(load_info), get_saved_data(load_info), path)
    end
    return configs_loaded
end

function load_configs(path) 
    configs = JSON3.read(path)
    configs = convert(Dict{Symbol, Any}, configs)
    load_dic_configs(configs, path)
end

function load_system_configs(path; ignore_info=false) 
    configs = JSON3.read(path)
    configs = convert(Dict{Symbol, Any}, configs)
    if ignore_info
        delete!(configs, :info)
    end
    sys_type = eval_saved_expression(configs[:sys_type])
    configs = load_dic_configs(configs, path)
    configs[:sys_type] = sys_type


    return configs
end

function load_state(path)
    load_info = get_load_info_serial(path)
    load_component(get_saved_type(load_info), get_saved_data(load_info), path)
end

function load_system(configs_loaded, rng, sys_type)
    System(
        state=configs_loaded[:state],
        space_cfg=configs_loaded[:space_cfg],
        dynamic_cfg=configs_loaded[:dynamic_cfg],
        int_cfg=configs_loaded[:int_cfg],
        info=configs_loaded[:info],
        time_info=configs_loaded[:time_info],
        debug_info=configs_loaded[:debug_info],
        rng=rng,
        sys_type=sys_type
    )
end

function load_system(configs_path::String, state_path::String, time_info_path=nothing)
    configs = JSON3.read(configs_path)
    configs = convert(Dict{Symbol, Any}, configs)
    configs[:state] = get_load_info_serial(state_path)
    sys_type = eval_saved_expression(configs[:sys_type])
    # sys_type = deserialize(joinpath(root, "sys_type.bin"))
    
    if !isnothing(time_info_path)
        configs[:time_info] = get_load_info_serial(time_info_path)
    end
    
    rng_file = joinpath(state_path, "rng.bin")
    if isfile(rng_file)
        rng = deserialize(joinpath(state_path, "rng.bin"))
    else
        rng = Random.GLOBAL_RNG
    end
    # rng = deserialize(joinpath(state_path, "rng.bin"))

    configs = load_dic_configs(configs, configs_path)
    load_system(configs, rng, sys_type)
end

load_system(root) = load_system(joinpath(root, "configs.json"), joinpath(root, "state"))

end # MaviSerder
