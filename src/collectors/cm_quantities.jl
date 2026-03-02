export CmQuantitiesCfg, CmsCfg, CmVelCfg

using Serialization, StaticArrays, StructTypes

using Mavi.Systems
using Mavi.States
using Mavi.Configs
using Mavi.MaviSerder
import Mavi.Experiments as Exp

abstract type CmQuantityCfg end
StructTypes.StructType(::Type{T}) where T <: CmQuantityCfg = StructTypes.Struct()

struct BaseCmQuantityState{T, N}
    data::Array{T, N}
end

function get_cm_quantity_shape(quantity::CmQuantityCfg) end
function get_cm_quantity_name(quantity::CmQuantityCfg) end

function get_base_cm_quantity_state(quantity::CmQuantityCfg, max_num_entities, num_collects, quantity_shape, NumT, system) 
    N = 2 + length(quantity_shape)
    data = Array{NumT, N}(undef, quantity_shape..., max_num_entities, num_collects)
    return BaseCmQuantityState(data)
end
get_cm_quantity_state(quantity::CmQuantityCfg, args...; kwargs...) = get_base_cm_quantity_state(quantity, args..., kwargs...)

function collect_cm_quantity(quantity::CmQuantityCfg, state, root_state, entities_ids, to_collect, system) end

get_save_data_cm_quantity(quantity::CmQuantityCfg, state) = state.data


struct CmsCfg <: CmQuantityCfg end
get_cm_quantity_shape(q::CmsCfg) = 2
get_cm_quantity_name(q::CmsCfg) = :cms


struct CmVelCfg <: CmQuantityCfg 
    num_dt_between::Int
end
get_cm_quantity_shape(q::CmVelCfg) = 2
get_cm_quantity_name(q::CmVelCfg) = :vel

mutable struct CmVelState{T, N} 
    cms2::Array{T, N}
    entities_ids::Vector{Int}
    num_frames::Int
    collect_cms2::Bool
end

function get_cm_quantity_state(quantity::CmVelCfg, args...; kwargs...)
    cms2 = get_base_cm_quantity_state(quantity, args..., kwargs...)    
    return CmVelState(cms2, Int[], 0, false)
end

function get_save_data_cm_quantity(quantity::CmVelCfg, state::CmVelState)
    return (cms2=state.cms2, num_frames=state.num_frames)
end

function collect_cm_quantity(quantity::CmVelCfg, state::CmVelState, root_state, entities_ids, to_collect, system) 
    if to_collect
        state.collect_cms2 = true
        state.entities_ids = entities_ids
    elseif state.collect_cms2 && system.time_info.num_steps - root_state.last_collect >= quantity.num_dt_between
        state.collect_cms2 = false
        data = state.cms2
        time_id = root_state.current_id
        for (idx, entity_id) in enumerate(state.entities_ids)
            data[:, idx, time_id] = get_cm(system, system.state, entity_id)
        end
        state.num_frames += 1
    end
end

@kwdef struct CmQuantitiesCfg{QC, M} <: Exp.ColCfg 
    delta_num_steps::Int
    max_num_entities::Int
    quantities_cfg::QC
    mask::M
end

function MaviSerder.get_obj_save_data(obj::CmQuantitiesCfg)
    if obj.mask isa Function
        obj = @set obj.mask = nothing
    end

    obj
end

mutable struct CmQuantitiesState{T, QS} <: Exp.ColState
    quantities_states::QS
    ids::Matrix{Int}
    num_entities::Vector{Int}
    times::Vector{T}
    last_collect::Int
    current_id::Int
    num_collects::Int
end

struct CmQuantitiesCol{T, QC, QS, M} <: Exp.Collector
    cfg::CmQuantitiesCfg{QC, M}
    state::CmQuantitiesState{T, QS}
end

function Exp.get_collector(col_cfg::CmQuantitiesCfg, exp_cfg::Exp.ExperimentCfg, system::System, state=nothing)
    if isnothing(state)
        q_states = []
        num_collects = floor(Int, (exp_cfg.tf - system.time_info.time) / system.int_cfg.dt / col_cfg.delta_num_steps)
        NumT = typeof(system.state).parameters[2]
        for q_cfg in col_cfg.quantities_cfg
            push!(q_states, get_cm_quantity_state(
                q_cfg, col_cfg.max_num_entities, num_collects, get_cm_quantity_shape(q_cfg), NumT, system
            ))
        end
        ids = Array{Int, 2}(undef, col_cfg.max_num_entities, num_collects)
        state = CmQuantitiesState(Tuple(q_states), ids, Int[], NumT[], system.time_info.num_steps, 1, num_collects)
    end

    CmQuantitiesCol(col_cfg, state)
end

get_mask_func(mask) = cm -> false
get_mask_func(mask::GeometryCfg) = cm -> !is_inside(cm, mask)
get_mask_func(mask::Function) = mask
    
function Exp.collect(col::CmQuantitiesCol, system)
    to_collect = false
    
    t_steps = system.time_info.num_steps
    state = col.state
    if !(t_steps - state.last_collect < col.cfg.delta_num_steps) && col.state.current_id <= col.state.num_collects
        to_collect = true
    end

    valid_entities_ids = Int[]
    num_valid = 0
    mask_func = get_mask_func(col.cfg.mask)
    if to_collect
        state.last_collect = t_steps

        cms_data = col.state.quantities_states[1].data
        current_time_id = col.state.current_id

        for entity_id in get_entities_ids(system)
            cm_i = get_cm(system, system.state, entity_id)
            if mask_func(cm_i)
                continue
            end
            
            push!(valid_entities_ids, entity_id)        
            num_valid += 1
            cms_data[:, num_valid, current_time_id] = cm_i
            
            if num_valid >= col.cfg.max_num_entities
                break
            end
        end
    end

    for (q, q_state) in zip(col.cfg.quantities_cfg[2:end], col.state.quantities_states[2:end])
        collect_cm_quantity(q, q_state, state, valid_entities_ids, to_collect, system)
    end

    if to_collect
        num_entities = length(valid_entities_ids)
        push!(state.times, system.time_info.time)
        push!(state.num_entities, num_entities)
        state.ids[1:num_entities, state.current_id] = valid_entities_ids 
        state.current_id += 1
    end
end

function Exp.save_data(col::CmQuantitiesCol, path)
    names = [get_cm_quantity_name(q) for q in col.cfg.quantities_cfg]
    datas = [get_save_data_cm_quantity(q_cfg, q_state) for (q_cfg, q_state) in zip(col.cfg.quantities_cfg, col.state.quantities_states)]
    all_datas = NamedTuple{Tuple(names)}(datas)

    final_data = (
        datas=all_datas,
        ids=col.state.ids,
        num_entities=col.state.num_entities,
        times=col.state.times,
        num_frames=col.state.current_id - 1,
    )
    serialize(joinpath(path, "data.bin"), final_data)
end

function Exp.load_data(::Type{CmQuantitiesCfg}, path)
    deserialize(joinpath(path, "data.bin"))
end
