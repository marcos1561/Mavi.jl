module Sources

export SourceCfg, Source
export SinkCfg, Sink
export update_area_empty!, process_sink_source!

using StaticArrays, StructTypes

import Mavi.Configs: GeometryCfg, RectangleCfg, check_intersection, is_inside
import Mavi.Rings: get_num_particles, get_ids

using Mavi.Rings.States
using Mavi.Rings.Configs
using Mavi.ChunksMod

struct RandomPol{T} end
StructTypes.StructType(::Type{RandomPol{T}}) where T = StructTypes.Struct()

get_spawn_pol(spawn_pol) = spawn_pol 
get_spawn_pol(spawn_pol::RandomPol{T}) where T = rand(T) * 2 * π 

struct SourceCfg{T<:AbstractFloat, SP}
    bottom_left::SVector{2, T}
    spawn_pos::Vector{SVector{2, T}}
    spawn_pol::SP
    pad::T
    offset::Tuple{T, T}
    size::Tuple{Int, Int}
end
function SourceCfg(; bottom_left, spawn_pos, spawn_pol, pad=0, offset=(0, 0), size=(1, 1))
    T = eltype(spawn_pos[1])
    offset_T = (convert(T, offset[1]), convert(T, offset[2]))
    
    if spawn_pol === :random
        spawn_pol = RandomPol{T}()
    end

    SourceCfg(bottom_left, spawn_pos, spawn_pol, pad, offset_T, size)
end

zero_offset(::Type{T}) where {T} = (zero(T), zero(T))

abstract type Checker end

struct ChunksChecker{P, C} <: Checker 
    pos::P
    chunks::C
    ids::Vector{Vector{CartesianIndex{2}}}
end
function ChunksChecker(pos, chunks, bbox_vec::Vector{G}) where G <: RectangleCfg 
    global_bbox = sum(bbox_vec)
    # global_bbox = RectangleCfg(
    #     length=global_bbox.length + cfg.pad*2,
    #     height=global_bbox.height + cfg.pad*2,
    #     bottom_left=global_bbox.bottom_left - SVector(cfg.pad, cfg.pad),
    # )
    bbox_ids = [CartesianIndex{2}[] for _ in 1:length(bbox_vec)]

    up = SVector(0, 1)
    chunk_tl = chunks.geometry_cfg.bottom_left + up * chunks.geometry_cfg.height

    found_start = false
    start_x_id, end_x_id = 0, 0
    for i in 1:chunks.num_cols
        x1 = chunk_tl[1] + (i-1) * chunks.chunk_length
        x2 = x1 + chunks.chunk_length
        if !found_start && x1 <= global_bbox.bottom_left[1] <= x2
            found_start = true
            start_x_id = i
        end

        if found_start && x1 >= global_bbox.bottom_left[1] + global_bbox.length
            end_x_id = i - 1
            break
        end
    end

    found_start = false
    start_y_id, end_y_id = 0, 0       
    for i in 1:chunks.num_rows
        y2 = chunk_tl[2] - (i-1) * chunks.chunk_height
        y1 = y2 - chunks.chunk_height
        if !found_start && y1 <= global_bbox.bottom_left[2] + global_bbox.height  <= y2
            found_start = true
            start_y_id = i
        end

        if found_start && y2 <= global_bbox.bottom_left[2]
            end_y_id = i - 1
            break
        end
    end

    for col in start_x_id:end_x_id
        for row in start_y_id:end_y_id
            rect = get_chunk_rect(chunks, row, col) 
            # RectangleCfg(
            #     length=chunks.chunk_length,
            #     height=chunks.chunk_height,
            #     bottom_left=chunk_tl - col*up + (row-1)*right,
            # )

            for (bbox, chunks_ids) in zip(bbox_vec, bbox_ids)
                if check_intersection(rect, bbox)
                    push!(chunks_ids, CartesianIndex(row, col))
                end
            end
        end
    end

    ChunksChecker(pos, chunks, bbox_ids)
end

struct PosChecker{P, PIDS} <: Checker
    pos::P
    part_ids::PIDS
end

struct Source{T<:AbstractFloat, C<:Checker, Cfg<:SourceCfg}
    cfg::Cfg
    checker::C
    is_empty::Vector{Bool}
    bbox::Vector{RectangleCfg{2, T}}
    bbox_spawn_pos::Vector{Vector{SVector{2, T}}}
end
function Source(cfg::SourceCfg{T, S}, chunks_checker_args=nothing, pos_checker_args=nothing) where {T, S}
    xs = [pos[1] for pos in cfg.spawn_pos]
    ys = [pos[2] for pos in cfg.spawn_pos]
    min_x = minimum(xs)
    max_x = maximum(xs)
    min_y = minimum(ys)
    max_y = maximum(ys)
    
    bbox_length = max_x - min_x + 2*cfg.pad
    bbox_height = max_y - min_y + 2*cfg.pad
    spawn_bbox = RectangleCfg(
        length=bbox_length, 
        height=bbox_height, 
        bottom_left=SVector(min_x, min_y),
    )

    bl = cfg.bottom_left
    right = SVector{2, T}(1, 0)
    up = SVector{2, T}(0, 1)
    bbox_vec::Vector{RectangleCfg{2, T}} = []
    bbox_pad_vec::Vector{RectangleCfg{2, T}} = []
    for i in 1:cfg.size[1]
        for j in 1:cfg.size[2]
            bl_i = bl + ((i-1)*bbox_length + i*cfg.offset[1]) * right  + ((j-1)*bbox_height + j*cfg.offset[2]) * up        
            bbox = RectangleCfg(
                length=bbox_length, 
                height=bbox_height, 
                bottom_left=bl_i,
            )
            bbox_pad = RectangleCfg(
                length=bbox_length+2*cfg.pad, 
                height=bbox_height+2*cfg.pad, 
                bottom_left=bl_i - SVector(cfg.pad, cfg.pad),
            )
            push!(bbox_vec, bbox)
            push!(bbox_pad_vec, bbox_pad)
        end
    end

    if !isnothing(chunks_checker_args)
        checker = ChunksChecker(chunks_checker_args..., bbox_pad_vec)
    elseif !isnothing(pos_checker_args)
        checker = PosChecker(pos_checker_args...)
    else
        error("No valid checker was provided! Please supply either `chunks` or `pos`.")
    end
    
    num_spawn = length(bbox_vec)
    spawn_pos = [SVector{2, T}[] for _ in 1:num_spawn]
    for idx in 1:num_spawn
        desloc = bbox_vec[idx].bottom_left - spawn_bbox.bottom_left + SVector(cfg.pad, cfg.pad)
        spawn_pos[idx] = cfg.spawn_pos .+ Scalar(desloc)
    end

    return Source(cfg, checker, fill(false, num_spawn), bbox_vec, spawn_pos)
end

function update_area_empty!(source::Source{T, C, Cfg}) where {T, C<:PosChecker, Cfg}
    checker = source.checker
    pos = checker.pos
    p_ids = get_ids(checker.part_ids)
    for idx in eachindex(source.bbox)
        source.is_empty[idx] = true
        bbox = source.bbox[idx]
        for p_idx in p_ids
            if is_inside(pos[p_idx], bbox, pad=source.cfg.pad)
                source.is_empty[idx] = false
                break
            end
        end
    end

    # for idx in eachindex(source.bbox)
    #     source.is_empty[idx] = true
    #     bbox = source.bbox[idx]
    #     found = false
    #     for ring_id in get_active_ids(checker.state)
    #         for p_id in 1:get_num_particles(dynamic_cfg, state, ring_id)
    #             p_scalar_idx = to_scalar_idx(state, ring_id, p_id)
    #             if is_inside(state.pos[p_scalar_idx], bbox, pad=source.cfg.pad)
    #                 source.is_empty[idx] = false
    #                 found = true
    #                 break
    #             end
    #         end
    #         if found
    #             break
    #         end
    #     end
    # end
end

function update_area_empty!(source::Source{T, C, Cfg}) where {T, C<:ChunksChecker, Cfg}
    checker = source.checker 
    chunks = checker.chunks
    pos = checker.pos
    for (i, chunks_ids) in enumerate(checker.ids)
        bbox = source.bbox[i]
        source.is_empty[i] = true
        for chunk_id in chunks_ids
            for p_idx in get_chunk_particles(chunks, chunk_id)
                if is_inside(pos[p_idx], bbox, pad=source.cfg.pad)
                    source.is_empty[i] = false
                    break
                end
            end
            if !source.is_empty[i]
                break
            end
        end
    end
end

function process_source!(source::Nothing, state=nothing) end

function process_source!(source, state, system=nothing)
    update_area_empty!(source)
    for idx in eachindex(source.is_empty)
        if !source.is_empty[idx]
            continue
        end
        ring_id = add_ring!(state, source.bbox_spawn_pos[idx], get_spawn_pol(source.cfg.spawn_pol))
        if !isnothing(system) && !isnothing(ring_id)
            num_p = get_num_particles(system, ring_id)
            system.info.cms[ring_id] = sum(state.rings_pos[1:num_p, ring_id]) / num_p
        end
    end
end

struct SinkCfg{G<:GeometryCfg}
    geometry_cfg::G
end

struct Sink{N, T, G}
    cfg::SinkCfg{G}
    pos::Vector{SVector{N, T}}
end

function process_sink!(sink::Sink, state)
    box = sink.cfg.geometry_cfg
    for ring_id in get_active_ids(state)
        if is_inside(sink.pos[ring_id], box)
            remove_ring!(state, ring_id)
        end
    end
end

process_sink_source!(sink_source::Sink, state, system=nothing) = process_sink!(sink_source, state)
process_sink_source!(sink_source::Source, state, system) = process_source!(sink_source, state, system)

end