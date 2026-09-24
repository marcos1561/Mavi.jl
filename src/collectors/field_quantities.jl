# ==
# Maskers
# ==
skip_entity(masker::Mavi.Configs.GeometryCfg, cm_i, entity_id, system) = !Mavi.Configs.is_inside(cm_i, masker)

# ==
# Field Quantity
# ==

abstract type FieldQuantityCfg end

function get_quantity_state(q_cfg, max_num_entities, num_collects, NumT, system)
    shape = get_quantity_shape(q_cfg, system)
    N = 2 + length(shape)
    return Array{NumT, N}(undef, shape..., max_num_entities, num_collects)
end

get_default_name(a::Type{FieldQuantityCfg}) = throw(ArgumentError("not implemented for $a"))

# function load_quantity_data(::Type{FieldQuantityCfg}, path)
#     return deserialize(joinpath(path, "data.bin"))
# end

struct PosCfg <: FieldQuantityCfg end
get_default_name(::Type{PosCfg}) = "pos"

function collect_quantity(q_cfg::PosCfg, q_state, state, valid_entities_ids, collect_flag, system)
    if !collect_flag
        return
    end

    current_id = state.current_id
    for (idx, i) in enumerate(valid_entities_ids)
        q_state[:, idx, current_id] .= Mavi.States.get_entity_pos(system, i)
    end
end

get_quantity_shape(q_cfg::PosCfg, system) = (length(system.state.pos[1]),)


struct VelCfg <: FieldQuantityCfg end
get_default_name(::Type{VelCfg}) = "vel"

function collect_quantity(q_cfg::VelCfg, q_state, state, valid_entities_ids, collect_flag, system)
    if !collect_flag
        return
    end

    current_id = state.current_id
    for (idx, i) in enumerate(valid_entities_ids)
        q_state[:, idx, current_id] .= Mavi.Systems.get_entity_vel(system, i)
    end
end

get_quantity_shape(q_cfg::VelCfg, system) = get_quantity_shape(PosCfg(), system)


struct PolCfg <: FieldQuantityCfg end
get_default_name(::Type{PolCfg}) = "pol"

function collect_quantity(q_cfg::PolCfg, q_state, state, valid_entities_ids, collect_flag, system)
    if !collect_flag
        return
    end

    current_id = state.current_id
    for (idx, i) in enumerate(valid_entities_ids)
        q_state[:, idx, current_id] .= Mavi.States.get_entity_pol(system.state, i)
    end
end

get_quantity_shape(q_cfg::PolCfg, system) = (1,)

# ==
# Field Quantities
# ==

struct FieldQuantitiesCfg{S, FG, M}  <: ColCfg
    sampler::S
    masker::M
    cfgs::FG
    names::Vector{String}
    max_num_entities::Int
end
function FieldQuantitiesCfg(; sampler, masker, cfgs, max_num_entities)
    if cfgs isa FieldQuantityCfg
        cfgs = Dict(get_default_name(typeof(cfgs))=>cfgs)
    elseif cfgs isa Vector
        cfgs = Dict(get_default_name(typeof(c)) => c for c in cfgs)
    end
    
    if !(cfgs isa Dict)
        throw(ArgumentError("`cfgs` must be one of the following: `FieldQuantityCfg`, `Vector{FieldQuantityCfg}` or  `Dict{String, FieldQuantityCfg}`, but it is $(typeof(cfgs))"))
    end

    @show cfgs

    cfg_list = Tuple(values(cfgs))
    names = string.(keys(cfgs))

    FieldQuantitiesCfg(sampler, masker, cfg_list, names, max_num_entities)
end

function get_collector(col_cfg::FieldQuantitiesCfg, exp_cfg::ExperimentCfg, system::System, state=nothing)
    if isnothing(state)
        q_states = []
        sampler = get_sampler(col_cfg.sampler)
        num_collects = get_num_collects(sampler.cfg, exp_cfg.tf)

        NumT = eltype(system.state.pos[1])
        for q_cfg in col_cfg.cfgs
            push!(q_states, get_quantity_state(
                q_cfg, col_cfg.max_num_entities, num_collects, NumT, system
            ))
        end
        ids = Array{Int, 2}(undef, col_cfg.max_num_entities, num_collects)
        
        state = FieldQuantitiesState(sampler, Tuple(q_states), ids, Int[], NumT[], 0, num_collects)
    end

    FieldQuantitiesCol(col_cfg, state)
end

mutable struct FieldQuantitiesState{S, QS, T}  <: ColState
    sampler::S
    quantities_states::QS
    ids::Matrix{Int}
    # valid_entities_ids::ValidEntitiesIds
    num_entities::Vector{Int}
    times::Vector{T}
    current_id::Int
    num_collects::Int
end

struct FieldQuantitiesCol{C, S} <: Collector
    cfg::C
    state::S
end

function collect(col::FieldQuantitiesCol, system; force_collect=false)
    time_info = system.time_info
    state = col.state

    collect_sampler = to_collect(state.sampler, time_info)
    collect_flag = (collect_sampler || force_collect) && (state.current_id < state.num_collects)

    if collect_sampler
        update_sampler(state.sampler, time_info)
    end

    if collect_flag
        state.current_id += 1
        num_valid = 0
        ids = @view state.ids[:, state.current_id]
        max_num_entities = col.cfg.max_num_entities
        for entity_id in Mavi.States.get_entities_ids(system)
            cm_i = Mavi.States.get_entity_pos(system, entity_id)
            if skip_entity(col.cfg.masker, cm_i, entity_id, system)
                continue
            end
            
            num_valid += 1
            ids[num_valid] = entity_id

            if num_valid >= max_num_entities
                break
            end
        end

        push!(state.times, time_info.time)
        push!(state.num_entities, num_valid)
    end

    current_id = state.current_id
    valid_entities_ids = @view state.ids[1:state.num_entities[current_id], current_id]
    for (q_cfg, q_state) in zip(col.cfg.cfgs, col.state.quantities_states)
        collect_quantity(q_cfg, q_state, state, valid_entities_ids, collect_flag, system)
    end

    # if collect_flag
    #     num_entities = length(valid_entities_ids)
    #     push!(state.times, system.time_info.time)
    #     push!(state.num_entities, num_entities)
    #     state.ids[1:num_entities, state.current_id] = valid_entities_ids 
    #     state.current_id += 1
    # end
end

function save_data(col::FieldQuantitiesCol, path)
    T_string = string(typeof(col.cfg))
    serialize(joinpath(path, "cfg_type.bin"), T_string)
    serialize(joinpath(path, "data.bin"), col)
end

function load_data(::Type{C}, path) where C <: FieldQuantitiesCfg
    deserialize(joinpath(path, "data.bin"))
end
