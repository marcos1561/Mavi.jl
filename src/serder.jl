module MaviSerder

export save_system, save_system_configs, save_component_serial, get_obj_save_data_json
export load_system, load_dic_configs, load_component, load_component_serial

using Serialization, JSON3
using StaticArrays

import Mavi
using Mavi.Systems

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
    data = get_obj_save_data_json(obj)
    serialize(joinpath(root_path, "type.bin"), data[:type])
    serialize(joinpath(root_path, "data.bin"), data[:data])
end

"Returns data to be saved inside `obj`."
get_obj_save_data(obj) = obj

"Returns data to be saved inside `obj` in the `json` format."
get_obj_save_data_json(obj) = (type=string(typeof(obj)), data=get_obj_save_data(obj))
get_obj_save_data_json(obj::Tuple) = [get_obj_save_data_json(x) for x in obj]

function save_system_configs(system::System, root)
    mkpath(root)
    configs_data = Dict(
        "sys_type" => system.type,
        "space_cfg" => get_obj_save_data_json(system.space_cfg),
        "dynamic_cfg" => get_obj_save_data_json(system.dynamic_cfg),
        "int_cfg" => get_obj_save_data_json(system.int_cfg),
        "info" => get_obj_save_data_json(system.info),
        "debug_info" => get_obj_save_data_json(system.debug_info),
        "time_info" => get_obj_save_data_json(system.time_info),
    )

    open(joinpath(root, "configs.json"), "w") do io
        JSON3.pretty(io, configs_data)
    end
end

"Save a system in the path `root`."
function save_system(system::System, root)
    save_system_configs(system, root)
    save_component_serial(system.state, root, "state")
end


# ==
# Deserialization
# ==

get_saved_type(load_info) = eval(Meta.parse(string(load_info[:type])))
get_saved_type(load_info::JSON3.Array) = Nothing

get_saved_data(load_info) = load_info[:data]
get_saved_data(load_info::JSON3.Array) = load_info

"Information to load an object saved using `Serialization`."
get_load_info_serial(root) = (type=deserialize(joinpath(root, "type.bin")), data=root)

"Load object of type `T` with data saved in `path`."
load_component(T::Type, path::String) = deserialize(joinpath(path, "data.bin"))

"Load object of type `T` with data in `data`."
load_component(T::Type, data::JSON3.Object) = JSON3.read(JSON3.write(data), T)

load_component(T::Type{Nothing}, _) = nothing

"Load array of objects with load information in `data`"
function load_component(T::Type{Nothing}, data::JSON3.Array)
    data_loaded = []
    for load_info in data
        push!(data_loaded, load_component(get_saved_type(load_info), get_saved_data(load_info)))
    end
    return Tuple(data_loaded)
end

function load_component_serial(root)
   load_info = get_load_info_serial(root)
   load_component(get_saved_type(load_info), get_saved_data(load_info)) 
end

"Load dictionary of object with load information in `configs`."
function load_dic_configs(configs)
    configs_loaded = Dict()
    for (name, load_info) in configs
        if name == :sys_type
            continue
        end
        configs_loaded[name] = load_component(get_saved_type(load_info), get_saved_data(load_info))
    end
    return configs_loaded
end

function load_system(configs, ::StandardSys)
    configs_loaded = load_dic_configs(configs)
    System(
        state=configs_loaded[:state],
        space_cfg=configs_loaded[:space_cfg],
        dynamic_cfg=configs_loaded[:dynamic_cfg],
        int_cfg=configs_loaded[:int_cfg],
        info=configs_loaded[:info],
        time_info=configs_loaded[:time_info],
        debug_info=configs_loaded[:debug_info],
    )
end

function load_system(configs_path::String, state_path::String, time_info_path=nothing)
    configs = JSON3.read(configs_path)
    configs = convert(Dict{Symbol, Any}, configs)
    configs[:state] = get_load_info_serial(state_path)
    sys_type = eval(Meta.parse(configs[:sys_type]))
    # sys_type = deserialize(joinpath(root, "sys_type.bin"))
    
    if !isnothing(time_info_path)
        configs[:time_info] = get_load_info_serial(time_info_path)
    end
    
    load_system(configs, sys_type)
end

load_system(root) = load_system(joinpath(root, "configs.json"), joinpath(root, "state"))

end # MaviSerder