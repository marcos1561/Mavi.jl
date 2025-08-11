module Rings

export RingsSystem, RingsState, ActiveState, Configs, InitStates, Integration
export get_num_total_particles, update_active_particles_ids!, update_part_ids!
export NeighborsCfg
export save, load_system

using StaticArrays, Serialization, JSON3, StructTypes

import Mavi
import Mavi.Systems: System, SystemType, get_num_total_particles, get_chunks, get_obj_save_data_json, load_dic_configs

include("states.jl")
include("configs.jl")
include("neighbors.jl")
using .States
using .Configs
using .NeighborsMod

struct RingsSys <: SystemType end

# 
# System Info
# 

@inline function Configs.get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:AbstractArray, T, I}
    return get_num_particles(dynamic_cfg, state.types[ring_id])
end

@inline function Configs.get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:Number, T, I}
    return dynamic_cfg.num_particles
end
@inline Configs.get_num_particles(system, ring_id) = Configs.get_num_particles(system.dynamic_cfg, system.state, ring_id)

@inline function Mavi.Systems.get_particle_radius(dynamic_cfg::RingsCfg{U, T, I}, state::RingsState, idx) where {U<:Number, T, I}
    return Configs.particle_radius(dynamic_cfg.interaction_finder)
end

@inline function Mavi.Systems.get_particle_radius(dynamic_cfg::RingsCfg{U, T, I}, state::RingsState, idx) where {U<:AbstractArray, T, I}
    type = state.types[get_ring_id(idx, num_max_particles(state))]
    inter = get_interaction_cfg(type, type, dynamic_cfg.interaction_finder)
    return Configs.particle_radius(inter)
end

function Mavi.Systems.particles_radius(dynamic_cfg::RingsCfg, state)
    p_radius::Vector{Float64} = []
    for ring_id in get_active_ids(state)
        inter = get_interaction_cfg(ring_id, ring_id, state, dynamic_cfg.interaction_finder)
        num_p = get_num_particles(dynamic_cfg, state, ring_id)
        radius_i = Mavi.Configs.particle_radius(inter)
        for _ in 1:num_p
            push!(p_radius, radius_i)
        end
    end

    return p_radius
end

function Mavi.Systems.is_valid_pair(state::RingsState, dynamic_cfg::RingsCfg, i, j)
    num_max = num_max_particles(state)
    for idx in (i, j)
        p_id, ring_id = get_particle_ring_id(idx, num_max)
        if !isnothing(state.active) && !state.active.mask[ring_id]
            return false
        end

        num_ring_particles = get_num_particles(dynamic_cfg, state, ring_id)
        if p_id > num_ring_particles
            return false
        end
    end
    return true
end

Mavi.Systems.get_particles_ids(system::System, state::RingsState) = get_ids(system.info.particles_ids) 

Mavi.Systems.get_num_total_particles(system::System, state::RingsState) = get_num(system.info.part_ids)

@kwdef mutable struct ParticleIds
    ids::Vector{Int}
    num::Int
end

get_ids(part_ids) = part_ids
function get_ids(part_ids::ParticleIds)
    ids = @view part_ids.ids[1:part_ids.num]
    return ids
end

get_num(part_ids) = length(part_ids)
get_num(part_ids::ParticleIds) = part_ids.num

function update_part_ids!(part_ids, state, dynamic_cfg) end
function update_part_ids!(part_ids::ParticleIds, state, dynamic_cfg)
    count = 1
    for ring_id in get_active_ids(state)
        first_id = to_scalar_idx(state, ring_id, 1)
        for i in 1:get_num_particles(dynamic_cfg, state, ring_id)
            part_ids.ids[count] = first_id +  i - 1
            count += 1
        end
    end
    part_ids.num = count - 1
end

function get_particles_ids_obj(state)
    if isnothing(state.types) && isnothing(state.active)
        return eachindex(state.pos)
    end

   ParticleIds(collect(1:length(state.pos)), 0)
end

include("sources.jl")
using .Sources


struct RingsInfo{T, PIDS, PN<:Union{ParticleNeighbors, Nothing}, RN<:Union{Neighbors, Nothing}, S, U}
    continuos_pos::Matrix{SVector{2, T}}
    areas::Vector{T}
    cms::Vector{SVector{2, T}}
    particles_ids::PIDS
    p_neigh::PN
    r_neigh::RN
    sources::S
    user_data::U
end
function RingsInfo(; state::RingsState, int_cfg, parts_ids,
    chunks=nothing, r_neighbors_cfg=nothing, p_neighbors_cfg=nothing, source_cfg=nothing, user_data=nothing)
    r_neighbors = nothing
    if !isnothing(r_neighbors_cfg) 
        num_max_neighbors = r_neighbors_cfg.only_count == false ? 15 : nothing 

        r_neighbors = Neighbors(
            num_entities=size(state.ring_pos, 3),
            num_max_neighbors=num_max_neighbors,
            device=int_cfg.device,
            cfg=r_neighbors_cfg,
        )   
    end

    p_neighbors = nothing
    if !isnothing(p_neighbors_cfg) 
        num_max_neighbors = p_neighbors_cfg.only_count == false ? 15 : nothing 

        neigh = Neighbors(
            num_entities=length(state.pos),
            num_max_neighbors=num_max_neighbors,
            device=int_cfg.device,
            cfg=p_neighbors_cfg,
        )   
        p_neighbors = ParticleNeighbors(
            neighbors=neigh,
            type=p_neighbors_cfg.type,
            num_max_particles=num_max_particles(state),
        )
    end
    
    cms = [zero(state.pos[1]) for _ in 1:size(state.rings_pos, 2)]
    
    sources = nothing
    if !isnothing(source_cfg)
        if !(hasmethod(iterate, (typeof(source_cfg),)))
            source_cfg = [source_cfg]
        end

        sources = [] 
        for s_cfg in source_cfg
            s = nothing
            if s_cfg isa SourceCfg
                if !isnothing(chunks)
                    s = Source(s_cfg, (state.pos, chunks), nothing)
                else
                    s = Source(s_cfg, nothing, (state.pos, parts_ids))
                end
            elseif s_cfg isa SinkCfg
                s = Sink(s_cfg, cms)
            else
                error("Unknown source configuration type: $(typeof(s_cfg))")
            end

            push!(sources, s)
        end
        sources = Tuple(sources)
    end

    RingsInfo(
        similar(state.rings_pos),
        Vector{Float64}(undef, size(state.rings_pos, 2)),
        cms,
        parts_ids,
        p_neighbors,
        r_neighbors,
        sources,
        user_data,
    )
end

Mavi.Systems.get_particles_ids(state::RingsState, info::RingsInfo) = get_ids(info.particles_ids)

function RingsSystem(;state, space_cfg, dynamic_cfg, int_cfg, p_neighbors_cfg=nothing,
    r_neighbors_cfg=nothing, source_cfg=nothing, user_data=nothing, time_info=nothing)
    if has_types_cfg(dynamic_cfg) != has_types_func(state)
        if has_types_cfg(dynamic_cfg)
            error("DynamicCfg has multiple types, but state.types is nothing!")
        else
            error("DynamicCfg has only one type, but state.types is not nothing!")
        end
    end

    # if !isnothing(r_neighbors) 
    #     num_max_neighbors = r_neighbors.only_count == false ? 15 : nothing 

    #     r_neighbors = Neighbors(
    #         num_entities=size(state.ring_pos, 3),
    #         num_max_neighbors=num_max_neighbors,
    #         device=int_cfg.device,
    #     )   
    # end

    # p_neighbors = nothing
    # if !isnothing(p_neighbors_cfg) 
    #     num_max_neighbors = p_neighbors_cfg.only_count == false ? 15 : nothing 

    #     neigh = Neighbors(
    #         num_entities=length(state.pos),
    #         num_max_neighbors=num_max_neighbors,
    #         device=int_cfg.device,
    #         cfg=p_neighbors_cfg,
    #     )   
    #     p_neighbors = ParticleNeighbors(
    #         neighbors=neigh,
    #         type=p_neighbors_cfg.type,
    #         num_max_particles=num_max_particles(state),
    #     )
    # end

    
    # parts_ids = get_particles_ids_obj(state)
    # update_part_ids!(parts_ids, state, dynamic_cfg)
    # cms = [zero(state.pos[1]) for _ in 1:size(state.rings_pos, 2)]
    
    # chunks = get_chunks(int_cfg, space_cfg, state, dynamic_cfg, parts_ids)

    # sources = nothing
    # if !isnothing(source_cfg)
    #     if !(typeof(source_cfg) <: AbstractArray)
    #         source_cfg = [source_cfg]
    #     end

    #     sources = [] 
    #     for s_cfg in source_cfg
    #         s = nothing
    #         if s_cfg isa SourceCfg
    #             if !isnothing(chunks)
    #                 s = Source(s_cfg, (state.pos, chunks), nothing)
    #             else
    #                 s = Source(s_cfg, nothing, (state.pos, parts_ids))
    #             end
    #         elseif s_cfg isa SinkCfg
    #             s = Sink(s_cfg, cms)
    #         else
    #             error("Unknown source configuration type: $(typeof(s_cfg))")
    #         end

    #         push!(sources, s)
    #     end
    #     sources = Tuple(sources)
    # end

    # info = RingsInfo(
    #     continuos_pos=similar(state.rings_pos),
    #     areas=Vector{Float64}(undef, size(state.rings_pos, 2)),
    #     cms=cms,
    #     particles_ids=parts_ids,
    #     p_neigh=p_neighbors, r_neigh=r_neighbors,
    #     sources=sources,
    #     user_data=user_data,
    # )

    parts_ids = get_particles_ids_obj(state)
    update_part_ids!(parts_ids, state, dynamic_cfg)
    chunks = get_chunks(int_cfg, space_cfg, state, dynamic_cfg, parts_ids)

    info = RingsInfo(
        state=state,
        int_cfg=int_cfg,
        parts_ids=parts_ids,
        chunks=chunks,
        r_neighbors_cfg=r_neighbors_cfg,
        p_neighbors_cfg=p_neighbors_cfg,
        source_cfg=source_cfg,
        user_data=user_data,
    )

    Mavi.System(
        state=state, 
        space_cfg=space_cfg, 
        dynamic_cfg=dynamic_cfg, 
        int_cfg=int_cfg,
        chunks=chunks, 
        info=info, 
        time_info=time_info,
        sys_type=RingsSys(),
    )
end

include("integration.jl")
include("utils.jl")
include("init_states.jl")

# include("view.jl")
# include("space_checks.jl")

# == 
# Serialization
# ==

using .Configs: InteractionCfg

StructTypes.StructType(::Type{Matrix{T}}) where T <: InteractionCfg = StructTypes.CustomStruct()

function StructTypes.lower(m::Matrix{T}) where T <: InteractionCfg
    return (size=size(m), data=vec(m))
end

function StructTypes.construct(::Type{Matrix{T}}, x) where T <: InteractionCfg
    elems = [isa(el, T) ? el : JSON3.read(JSON3.write(el), T) for el in x["data"]]
    reshape(elems, x["size"]...)
end

function Mavi.Systems.get_obj_save_data(info::RingsInfo)
    r_neigh_cfg = nothing
    if !isnothing(info.r_neigh)
        r_neigh_cfg = info.r_neigh.cfg
    end
    # serialize(joinpath(info_root, "r_neighbors_cfg.bin"), info.r_neigh.cfg)
    
    p_neigh_cfg = nothing
    if !isnothing(info.p_neigh)
        p_neigh_cfg = info.p_neigh.cfg
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

function Mavi.Systems.load_component(T::Type{RI}, data::JSON3.Object) where RI <: RingsInfo 
    # configs_loaded = Dict()
    # println("== Info ==")
    # for (name, value) in data
    #     println("===")
    #     println(name)
    #     println(value[:type])
    #     println(typeof(value[:data]))
        
    #     T = eval(Meta.parse(string(value[:type]))) 
    #     configs_loaded[name] = Mavi.Systems.load_component(T, value[:data])
    # end
    # return configs_loaded
    return load_dic_configs(data)
end

function Mavi.Systems.get_obj_save_data(state::RingsState)
    if isnothing(state.active)
        active_mask = nothing
    else
        active_mask = state.active.mask
    end
    return (
        pol=state.pol, 
        rings_pos=state.rings_pos,
        types=state.types,
        active_state=active_mask,
    )
end

# function Mavi.Systems.save_component_serial(state::RingsState, root, name)
#     root_path = mkpath(joinpath(root, name))
#     data = get_obj_type_and_data(obj)
#     serialize(joinpath(root_path, "type.bin"), data[:type])
#     serialize(joinpath(root_path, "data.bin"), data[:data])
# end

function Mavi.Systems.load_component(T::Type{RS}, path) where RS <: RingsState 
    data = deserialize(joinpath(path, "data.bin"))
    
    active_state = isnothing(data.active_state) ? nothing : ActiveState(data.active_state)

    RingsState(
        rings_pos=data.rings_pos,
        pol=data.pol,
        types=data.types,
        active_state=active_state,
    )
end

function Mavi.Systems.load_system(configs, ::RingsSys)
    configs_loaded = load_dic_configs(configs)
    # configs_loaded = Dict()
    # for (name, value) in configs
    #     if name == :sys_type
    #         continue
    #     end
        
    #     println("=====")
    #     println(name)
    #     println(typeof(value[:data]))

    #     T = eval(Meta.parse(string(value[:type]))) 
    #     configs_loaded[name] = Mavi.Systems.load_component(T, value[:data])
    # end

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
    )

end

# function Mavi.Systems.save_comp(info::RingsInfo, root, name)
#     info_root = joinpath(root, name)
#     mkpath(info_root)
    
#     r_neigh_cfg = nothing
#     if !isnothing(info.r_neigh)
#         r_neigh_cfg = info.r_neigh.cfg
#     end
#     serialize(joinpath(info_root, "r_neighbors_cfg.bin"), info.r_neigh.cfg)
    
#     p_neigh_cfg = nothing
#     if !isnothing(info.p_neigh)
#         p_neigh_cfg = info.p_neigh.cfg
#     end
#     serialize(joinpath(info_root, "p_neighbors_cfg.bin"), p_neigh_cfg)
    
#     source_cfg = nothing
#     if !isnothing(info.sources)
#         source_cfg = Tuple([s.cfg for s in info.sources])
#     end
#     serialize(joinpath(info_root, "source_cfg.bin"), source_cfg)     
    
#     Mavi.Systems.save(info.user_data, info_root, "user_data")
# end

# function Mavi.Systems.load_system(root, ::RingsSys)
#     main_members = load_main_members(root)

#     info_root = joinpath(root, "info") 

#     p_neighbors_cfg = deserialize(joinpath(info_root, "p_neighbors_cfg.bin"))
#     r_neighbors_cfg = deserialize(joinpath(info_root, "r_neighbors_cfg.bin"))
#     source_cfg = deserialize(joinpath(info_root, "source_cfg.bin"))
#     user_data = deserialize(joinpath(info_root, "user_data.bin"))
    
#     RingsSystem(
#         state=main_members.state,
#         space_cfg=main_members.space_cfg,
#         dynamic_cfg=main_members.dynamic_cfg,
#         int_cfg=main_members.int_cfg,
#         p_neighbors_cfg=p_neighbors_cfg,
#         r_neighbors_cfg=r_neighbors_cfg,
#         source_cfg=source_cfg,
#         user_data=user_data,
#         time_info=main_members.time_info,
#     )
# end

end