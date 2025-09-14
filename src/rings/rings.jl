module Rings

export RingsSystem, RingsState, ActiveState, Configs, InitStates, Integration
export get_num_particles, get_num_total_particles
export update_active_particles_ids!, update_part_ids!
export get_continuos_pos, ring_points
export NeighborsCfg, Invasion
export save_system, load_system

using StaticArrays

using Mavi
import Mavi.Systems: System, SystemType, get_num_total_particles, get_chunks
import Mavi.Configs: SpaceCfg, PeriodicWalls, ManyWalls
import Mavi.ChunksMod: Chunks, update_chunks!

include("states.jl")
include("configs.jl")
include("neighbors.jl")
using .States
using .Configs
using .NeighborsMod

struct RingsSys <: SystemType end

# =
# System utils functions 
# =

@inline function Configs.get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:AbstractArray, T, I}
    return get_num_particles(dynamic_cfg, state.types[ring_id])
end

@inline function Configs.get_num_particles(dynamic_cfg::RingsCfg{U, T, I}, state, ring_id) where {U<:Number, T, I}
    return dynamic_cfg.num_particles
end
@inline Configs.get_num_particles(system, ring_id) = Configs.get_num_particles(system.dynamic_cfg, system.state, ring_id)

function get_continuos_pos(ring_id, system, wall_type)
    num_particles = get_num_particles(system.dynamic_cfg, system.state, ring_id) 
    return @view system.state.rings_pos[1:num_particles, ring_id]
end

function get_continuos_pos(ring_id, system, wall_type::PeriodicWalls)
    num_particles = get_num_particles(system.dynamic_cfg, system.state, ring_id) 
    return @view system.info.continuos_pos[1:num_particles, ring_id]
end

get_continuos_pos(ring_id, system, wall_type::ManyWalls) = get_continuos_pos(ring_id, system, wall_type.list[1])

ring_points(system, ring_id) = get_continuos_pos(ring_id, system, system.space_cfg.wall_type)

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

# =
# Particles Indices
# = 

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

# =
# Invasions
# = 

@kwdef struct Invasion
    invasor::Int
    invaded::Int
    p_id::Int
end

mutable struct InvasionsInfo
    list::Vector{Invasion}
    last_check::Int
end
InvasionsInfo() = InvasionsInfo(Invasion[], 0)
InvasionsInfo(last_check::Int) = InvasionsInfo(Invasion[], last_check)

# =
# Rings Info
# =

include("sources.jl")
using .Sources

struct RingsInfo{T, PIDS, PN<:Union{ParticleNeighbors, Nothing}, RN<:Union{Neighbors, Nothing}, RC<:Union{Chunks, Nothing}, S, U}
    continuos_pos::Matrix{SVector{2, T}}
    areas::Vector{T}
    cms::Vector{SVector{2, T}}
    particles_ids::PIDS
    p_neigh::PN
    r_neigh::RN
    r_chunks::RC
    invasions::InvasionsInfo
    sources::S
    user_data::U
end
function RingsInfo(; dynamic_cfg, space_cfg, state::RingsState, int_cfg, parts_ids,
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

    rings_chunks = nothing
    if !isnothing(int_cfg.extra.r_chunks_cfg)
        r_chunks_cfg = int_cfg.extra.r_chunks_cfg

        bounding_box = Configs.get_bounding_box(space_cfg.geometry_cfg) 
        chunks_space_cfg = SpaceCfg(
            wall_type=space_cfg.wall_type,
            geometry_cfg=bounding_box,
        )

        rings_r = minimum(get_ring_radius(dynamic_cfg))

        rings_chunks = Chunks(
            r_chunks_cfg.num_cols, r_chunks_cfg.num_rows,
            chunks_space_cfg, cms, rings_r, extra_info=state,
        )
    end

    RingsInfo(
        similar(state.rings_pos),
        Vector{Float64}(undef, size(state.rings_pos, 2)),
        cms,
        parts_ids,
        p_neighbors,
        r_neighbors,
        rings_chunks,
        InvasionsInfo(),
        sources,
        user_data,
    )
end

Mavi.Systems.get_particles_ids(state::RingsState, info::RingsInfo) = get_ids(info.particles_ids)

include("integration.jl")
include("utils.jl")
include("init_states.jl")
include("serder.jl")
include("collectors.jl")

# =
# Rings System
# = 

function RingsSystem(;state, space_cfg, dynamic_cfg, int_cfg, p_neighbors_cfg=nothing,
    r_neighbors_cfg=nothing, source_cfg=nothing, user_data=nothing, time_info=nothing)
    if has_types_cfg(dynamic_cfg) != has_types_func(state)
        if has_types_cfg(dynamic_cfg)
            error("DynamicCfg has multiple types, but state.types is nothing!")
        else
            error("DynamicCfg has only one type, but state.types is not nothing!")
        end
    end

    parts_ids = get_particles_ids_obj(state)
    update_part_ids!(parts_ids, state, dynamic_cfg)
    chunks = get_chunks(int_cfg.chunks_cfg, space_cfg, state.pos, dynamic_cfg, parts_ids)

    info = RingsInfo(
        dynamic_cfg=dynamic_cfg,
        space_cfg=space_cfg,
        state=state,
        int_cfg=int_cfg,
        parts_ids=parts_ids,
        chunks=chunks,
        r_neighbors_cfg=r_neighbors_cfg,
        p_neighbors_cfg=p_neighbors_cfg,
        source_cfg=source_cfg,
        user_data=user_data,
    )

    system = Mavi.System(
        state=state, 
        space_cfg=space_cfg, 
        dynamic_cfg=dynamic_cfg, 
        int_cfg=int_cfg,
        chunks=chunks, 
        info=info, 
        time_info=time_info,
        sys_type=RingsSys(),
    )

    if !isnothing(system.int_cfg.extra.invasions_cfg)
        info.invasions.last_check = system.time_info.num_steps
    end

    Integration.update_continuos_pos!(system, system.space_cfg.wall_type)
    Integration.update_ids!(system)
    Integration.update_cms!(system)
    Integration.update_chunks_all!(system)
    Integration.update_invasions!(system) 
    Integration.check_invasions!(system)

    Integration.cleaning!(system)
    Integration.forces!(system)
    Integration.neigh_sum_buffers(system.info.p_neigh)

    return system
end

end