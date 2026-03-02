module RingsCollectors

export RingQuantitiesCfg, RingCmsCfg, RingVelCfg

using Serialization, StaticArrays

using Mavi.Systems
using Mavi.Rings.States
import Mavi.Experiments as Exp

abstract type RingQuantityCfg end

struct BaseRingQuantityState{T, N}
    data::Array{T, N}
end

function get_ring_quantity_shape(quantity::RingQuantityCfg) end
function get_ring_quantity_name(quantity::RingQuantityCfg) end

function get_base_ring_quantity_state(quantity::RingQuantityCfg, max_num_rings, num_collects, quantity_shape, NumT, system) 
    data = Array{NumT, 2 + length(quantity_shape)}(undef, quantity_shape..., max_num_rings, num_collects)
    return BaseRingQuantityState(data)
end
get_ring_quantity_state(quantity::RingQuantityCfg, args...; kwargs...) = get_base_ring_quantity_state(quantity, args..., kwargs...)

function collect_ring_quantity(quantity::RingQuantityCfg, state, root_state, rings_ids, to_collect, system) end

get_save_data_ring_quantity(quantity::RingQuantityCfg, state) = state.data


struct RingCmsCfg <: RingQuantityCfg end
get_ring_quantity_shape(q::RingCmsCfg) = 2
get_ring_quantity_name(q::RingCmsCfg) = :cms


struct RingVelCfg <: RingQuantityCfg 
    num_dt_between::Int
end
get_ring_quantity_shape(q::RingVelCfg) = 2
get_ring_quantity_name(q::RingVelCfg) = :vel

mutable struct RingsVelState{T, N} 
    cms2::Array{T, N}
    rings_ids::Vector{Int}
    num_frames::Int
    collect_cms2::Bool
end

function get_ring_quantity_state(quantity::RingVelCfg, args...; kwargs...)
    cms2 = get_base_ring_quantity_state(quantity, max_num_rings, args..., kwargs...)    
    return RingsVelState(cms2, Int[], 0, false)
end

function get_save_data_ring_quantity(quantity::RingVelCfg, state::RingsVelState)
    return (cms2=state.cms2, num_frames=state.num_frames)
end

function collect_ring_quantity(quantity::RingVelCfg, state::RingsVelState, root_state, rings_ids, to_collect, system) 
    if to_collect
        state.collect_cms2 = true
        state.rings_ids = rings_ids
    elseif state.collect_cms2 && system.time_info.num_steps - root_state.last_collect >= quantity.num_dt_between
        state.collect_cms2 = false
        data = state.cms2
        time_id = root_state.current_id
        for (idx, ring_id) in enumerate(state.rings_ids)
            data[:, idx, time_id] = system.info.cms[ring_id]
        end
        state.num_frames += 1
    end
end


@kwdef struct RingQuantitiesCfg{QC, F} <: Exp.ColCfg 
    delta_num_steps::Int
    max_num_rings::Int
    quantities_cfg::QC
    mask::F
end

mutable struct RingQuantitiesState{T, QS} <: Exp.ColState
    quantities_states::QS
    ids::Matrix{Int}
    num_rings::Vector{Int}
    times::Vector{T}
    last_collect::Float64
    current_id::Int
    num_collects::Int
end

struct RingQuantitiesCol{T, QC, QS} <: Exp.Collector
    cfg::RingQuantitiesCfg{QC}
    state::RingQuantitiesState{T, QS}
end

function Exp.get_collector(col_cfg::RingQuantitiesCfg, exp_cfg::Exp.ExperimentCfg, system::System, state=nothing)
    if isnothing(state)
        q_states = []
        num_collects =  floor(Int, (exp_cfg.tf - system.time_info.time) / system.int_cfg.dt)
        for q_cfg in col_cfg.quantities_cfg
            NumT = typeof(system.state).parameters[1]
            push!(q_states, get_ring_quantity_state(
                q_cfg, col_cfg.max_num_rings, num_collects, get_ring_quantity_shape(q_cfg), NumT, system
            ))
        end
        ids = Array{Int, 2}(undef, col_cfg.max_num_rings, num_collects)
        state = RingQuantitiesState(Tuple(q_states...), ids, Int[], NumT[], system.time_info.num_steps, 1, num_collects)
    end
    RingQuantitiesCol(col_cfg, state)
end

function Exp.collect(col::RingQuantitiesCol, system)
    to_collect = false
    
    t_steps = system.time_info.num_steps
    state = col.state
    if !(t_steps - state.last_collect < col.cfg.num_steps) && col.state.current_id <= col.state.num_collects
        to_collect = true
    end

    valid_rings_ids = Int[]
    num_valid = 0
    if to_collect
        state.last_collect = t_steps

        cms_data = col.state.quantities_states[1]
        current_time_id = col.state.current_id

        for ring_id in get_rings_ids(system)
            cm_i = system.info.cms[ring_id]
            if col.cfg.mask(cm_i)
                continue
            end
            
            push!(valid_rings_ids, ring_id)        
            num_valid += 1
            cms_data[:, num_valid, current_time_id] = cm_i
            
            if num_valid >= col.cfg.max_num_rings
                break
            end
        end

        push!(state.num_rings, length(valid_rings_ids))
    end

    for (q, q_state) in zip(col.cfg.quantities_cfg[2:end], col.state.quantities_states[2:end])
        collect_ring_quantity(q, q_state, state, valid_rings_ids, to_collect, system)
    end

    if to_collect
        col.state.current_id += 1
    end
end

function Exp.save_data(col::RingQuantitiesCol, path)
    names = [get_ring_quantity_name(q) for q in col.cfg.quantities_cfg]
    datas = [get_save_data_ring_quantity(q_cfg, q_state) for (q_cfg, q_state) in zip(col.cfg.quantities_cfg, col.state.quantities_states)]
    all_datas = NamedTuple{Tuple(names)}(datas)

    final_data = (
        datas=all_datas,
        ids=col.state.ids,
        num_rings=col.state.num_rings,
        times=col.state.times,
        num_frames=col.state.current_id - 1,
    )
    serialize(joinpath(path, "data.bin"), final_data)
end

function Exp.load_data(::Type{RingQuantitiesCfg}, path)
    deserialize(joinpath(path, "data.bin"))
end


end # RingsCollectors