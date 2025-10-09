module Rings

export RingsSystem, RingsState, RingsSys, ActiveState, Configs, InitStates, Integration
export get_num_particles, get_num_total_particles
export get_continuos_pos, ring_points
export NeighborsCfg, Invasion
export save_system, load_system

using StaticArrays

using Mavi
import Mavi.Systems: System, SystemType, get_num_total_particles, get_chunks, system_deepcopy
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

States.ring_num_particles(system, ring_id) = ring_num_particles(system.state, ring_id)

function get_continuos_pos(system, ring_id, wall_type)
    num_particles = ring_num_particles(system.state, ring_id) 
    return @view system.state.rings_pos[1:num_particles, ring_id]
end

function get_continuos_pos(system, ring_id, wall_type::PeriodicWalls)
    num_particles = ring_num_particles(system.state, ring_id) 
    return @view system.info.continuos_pos[1:num_particles, ring_id]
end

get_continuos_pos(system, ring_id, wall_type::ManyWalls) = get_continuos_pos(system, ring_id, wall_type.list[1])

ring_points(system, ring_id) = get_continuos_pos(system, ring_id, system.space_cfg.wall_type)

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
    for ring_id in get_rings_ids(state)
        inter = get_interaction_cfg(ring_id, ring_id, state, dynamic_cfg.interaction_finder)
        num_p = ring_num_particles(state, ring_id)
        radius_i = Mavi.Configs.particle_radius(inter)
        for _ in 1:num_p
            push!(p_radius, radius_i)
        end
    end

    return p_radius
end

function Mavi.Systems.is_valid_pair(state::RingsState, i, j)
    num_max = num_max_particles(state)
    for idx in (i, j)
        p_id, ring_id = get_particle_ring_id(idx, num_max)
        if !isnothing(state.active) && !state.active.mask[ring_id]
            return false
        end

        num_ring_particles = ring_num_particles(state, ring_id)
        if p_id > num_ring_particles
            return false
        end
    end
    return true
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
# Chunk Info
# = 

struct RingsChunksInfo{S, F}
    state::S
    ids_func::F
end

# =
# Rings Info
# =

include("sources.jl")
using .Sources

struct RingsInfo{T, PN<:Union{ParticleNeighbors, Nothing}, RN<:Union{Neighbors, Nothing}, RC<:Union{Chunks, Nothing}, S, U}
    continuos_pos::Matrix{SVector{2, T}}
    areas::Vector{T}
    cms::Vector{SVector{2, T}}
    p_neigh::PN
    r_neigh::RN
    r_chunks::RC
    invasions::InvasionsInfo
    sources::S
    user_data::U
end
function RingsInfo(; dynamic_cfg, space_cfg, state::RingsState, int_cfg,
    chunks=nothing, r_neighbors_cfg=nothing, p_neighbors_cfg=nothing, source_cfg=nothing, user_data=nothing)
    r_neighbors = nothing
    if !isnothing(r_neighbors_cfg) 
        num_max_neighbors = r_neighbors_cfg.only_count == false ? 15 : nothing 

        r_neighbors = Neighbors(
            num_entities=size(state.rings_pos, 2),
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
                    s = Source(s_cfg, nothing, (state.pos, state.rings_ids))
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
            chunks_space_cfg, cms, rings_r, extra_info=RingsChunksInfo(state, States.get_rings_ids),
        )
    end

    RingsInfo(
        similar(state.rings_pos),
        Vector{Float64}(undef, size(state.rings_pos, 2)),
        cms,
        p_neighbors,
        r_neighbors,
        rings_chunks,
        InvasionsInfo(),
        sources,
        user_data,
    )
end

# Mavi.Systems.get_particles_ids(state::RingsState, info::RingsInfo) = get_ids(info.particles_ids)

include("integration.jl")
include("utils.jl")
include("init_states.jl")
include("serder.jl")
include("collectors.jl")

# =
# Rings System
# = 

function RingsSystem(;state, space_cfg, dynamic_cfg, int_cfg, p_neighbors_cfg=nothing,
    r_neighbors_cfg=nothing, source_cfg=nothing, user_data=nothing, time_info=nothing, rng=nothing)
    if has_types_cfg(dynamic_cfg) != has_types_func(state)
        if has_types_cfg(dynamic_cfg)
            error("DynamicCfg has multiple types, but state.types is nothing!")
        else
            error("DynamicCfg has only one type, but state.types is not nothing!")
        end
    end

    chunks = get_chunks(int_cfg.chunks_cfg, space_cfg, state.pos, dynamic_cfg, RingsChunksInfo(state, States.get_particles_ids))
    # chunks = get_chunks(int_cfg.chunks_cfg, space_cfg, state.pos, dynamic_cfg, state)

    info = RingsInfo(
        dynamic_cfg=dynamic_cfg,
        space_cfg=space_cfg,
        state=state,
        int_cfg=int_cfg,
        chunks=chunks,
        r_neighbors_cfg=r_neighbors_cfg,
        p_neighbors_cfg=p_neighbors_cfg,
        source_cfg=source_cfg,
        user_data=user_data,
    )
    
    if dynamic_cfg.num_particles == -1
        dynamic_cfg = RingsCfg(dynamic_cfg, state.num_particles)
    end

    if dynamic_cfg.num_particles != state.num_particles
        error("`dynamic_cfg.num_particles` is not equal to `state.num_particles`!")
    end

    system = Mavi.System(
        state=state, 
        space_cfg=space_cfg, 
        dynamic_cfg=dynamic_cfg, 
        int_cfg=int_cfg,
        chunks=chunks, 
        info=info, 
        time_info=time_info,
        sys_type=RingsSys(),
        rng=rng,
    )

    if !isnothing(system.int_cfg.extra.invasions_cfg)
        info.invasions.last_check = system.time_info.num_steps
    end

    Integration.update_continuos_pos!(system, system.space_cfg.wall_type)
    Integration.update_ids!(system)
    Integration.update_cms!(system)
    Integration.update_chunks_all!(system)
    Integration.update_invasions!(system) 

    Integration.cleaning!(system)
    Integration.forces!(system)
    Integration.neigh_sum_buffers(system.info.p_neigh)

    return system
end

function system_deepcopy(system::System, sys_type::RingsSys)
    sources = system.info.sources

    source_cfg = isnothing(sources) ? nothing : [deepcopy(s.cfg) for s in system.info.sources]
    
    RingsSystem(
        state=deepcopy(system.state),
        space_cfg=deepcopy(system.space_cfg),
        dynamic_cfg=deepcopy(system.dynamic_cfg),
        int_cfg=deepcopy(system.int_cfg),
        p_neighbors_cfg=isnothing(system.info.p_neigh) ? nothing : deepcopy(get_neigh(system.info.p_neigh).cfg),
        r_neighbors_cfg=isnothing(system.info.r_neigh) ? nothing : deepcopy(get_neigh(system.info.r_neigh).cfg),
        source_cfg=source_cfg,
        user_data=deepcopy(system.info.user_data),
        time_info=deepcopy(system.time_info),
        rng=deepcopy(system.rng),
    )
end

end