module RingsSerder

using Serialization, JSON3, StructTypes

import Mavi.MaviSerder: MaviSerder, get_obj_save_data_json, load_dic_configs

using Mavi.Rings
using Mavi.Rings: RingsInfo, RingsSys
using Mavi.Rings.Configs: InteractionCfg

# ==
# Matrix of InteractionCfg
# ==

StructTypes.StructType(::Type{Matrix{T}}) where T <: InteractionCfg = StructTypes.CustomStruct()

function StructTypes.lower(m::Matrix{T}) where T <: InteractionCfg
    return (size=size(m), data=vec(m))
end

function StructTypes.construct(::Type{Matrix{T}}, x) where T <: InteractionCfg
    elems = [isa(el, T) ? el : JSON3.read(JSON3.write(el), T) for el in x["data"]]
    reshape(elems, x["size"]...)
end

# ==
# RingsInfo
# ==

function MaviSerder.get_obj_save_data(info::RingsInfo)
    r_neigh_cfg = nothing
    if !isnothing(info.r_neigh)
        r_neigh_cfg = info.r_neigh.cfg
    end
    # serialize(joinpath(info_root, "r_neighbors_cfg.bin"), info.r_neigh.cfg)
    
    p_neigh_cfg = nothing
    if !isnothing(info.p_neigh)
        p_neigh_cfg = info.p_neigh.neighbors.cfg
    end
    # serialize(joinpath(info_root, "p_neighbors_cfg.bin"), p_neigh_cfg)
    
    source_cfg = nothing
    if !isnothing(info.sources)
        source_cfg = Tuple([s.cfg for s in info.sources])
    end
    # serialize(joinpath(info_root, "source_cfg.bin"), source_cfg)     
    
    # Mavi.Systems.save(info.user_data, info_root, "user_data")

    return (
        p_neighbors_cfg=get_obj_save_data_json(p_neigh_cfg), 
        r_neighbors_cfg=get_obj_save_data_json(r_neigh_cfg), 
        source_cfg=get_obj_save_data_json(source_cfg), 
        user_data=get_obj_save_data_json(info.user_data),
    )
end

function MaviSerder.load_component(T::Type{RI}, data::JSON3.Object) where RI <: RingsInfo 
    return load_dic_configs(data)
end

# ==
# RingsState
# ==

function MaviSerder.get_obj_save_data(state::RingsState)
    if state.rings_ids isa Rings.States.VarRingsIds
        active_mask = state.rings_ids.mask
    else
        active_mask = nothing
    end
    return (
        pol=state.pol, 
        rings_pos=state.rings_pos,
        types=state.types,
        num_particles=state.num_particles,
        active_state=active_mask,
    )
end

function MaviSerder.load_component(T::Type{RS}, path::String) where RS <: RingsState 
    data = deserialize(joinpath(path, "data.bin"))
    
    active_state = isnothing(data.active_state) ? nothing : ActiveState(data.active_state)

    RingsState(
        rings_pos=data.rings_pos,
        pol=data.pol,
        types=data.types,
        active_state=active_state,
        num_particles=data.num_particles,
    )
end

# ==
# Rings System
# ==

function MaviSerder.load_system(configs, rng, ::RingsSys)
    configs_loaded = load_dic_configs(configs)

    RingsSystem(
        state=configs_loaded[:state],
        space_cfg=configs_loaded[:space_cfg],
        dynamic_cfg=configs_loaded[:dynamic_cfg],
        int_cfg=configs_loaded[:int_cfg],
        p_neighbors_cfg=configs_loaded[:info][:p_neighbors_cfg],
        r_neighbors_cfg=configs_loaded[:info][:r_neighbors_cfg],
        source_cfg=configs_loaded[:info][:source_cfg],
        user_data=configs_loaded[:info][:user_data],
        time_info=configs_loaded[:time_info],
        rng=rng,
    )

end


end # RingsSerder