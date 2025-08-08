module Sources

export SourceCfg, Source
export update_area_empty!, process_source!

using StaticArrays

import Mavi.Configs: RectangleCfg, check_intersection, is_inside
import Mavi.Rings: get_num_particles, get_ids

using Mavi.Rings.States
using Mavi.Rings.Configs
using Mavi.ChunksMod

struct RandomPol{T} end

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

struct ChunksChecker{C} <: Checker 
    chunks::C
    ids::Vector{Vector{CartesianIndex{2}}}
end

struct PosChecker{P, PIDS} <: Checker
    pos::P
    part_ids::PIDS
end

struct Source{T<:AbstractFloat, C<:Checker}
    cfg::SourceCfg
    checker::C
    is_empty::Vector{Bool}
    bbox::Vector{RectangleCfg{2, T}}
    bbox_spawn_pos::Vector{Vector{SVector{2, T}}}
end
function Source(cfg::SourceCfg{T, S}, chunks::Union{Chunks, Nothing}=nothing, pos_checker_args=nothing) where {T, S}
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
    for i in 1:cfg.size[1]
        for j in 1:cfg.size[2]
            bl_i = bl + ((i-1)*bbox_length + i*cfg.offset[1]) * right  + ((j-1)*bbox_height + j*cfg.offset[2]) * up        
            bbox = RectangleCfg(
                length=bbox_length, 
                height=bbox_height, 
                bottom_left=bl_i,
            )
            push!(bbox_vec, bbox)
        end
    end

    if !isnothing(chunks)
        global_bbox = sum(bbox_vec)
        global_bbox = RectangleCfg(
            length=global_bbox.length + cfg.pad*2,
            height=global_bbox.height + cfg.pad*2,
            bottom_left=global_bbox.bottom_left - SVector(cfg.pad, cfg.pad),
        )
        bbox_ids = [CartesianIndex{2}[] for _ in 1:length(bbox_vec)]

        chunk_tl = chunks.geometry_cfg.bottom_left + up * chunks.geometry_cfg.height

        found_start = false
        for i in 1:chunks.num_cols
            x1 = chunk_tl[1] + (i-1) * chunks.chunk_length
            x2 = x1 + chunks.chunk_length
            if !found_start && x1 <= global_bbox.bottom_left[1] <= x2
                found_start = true
                start_x_id = i
            end

            if found_start && x1 > global_bbox.bottom_left[1] + global_bbox.length
                end_x_id = i - 1
                break
            end
        end

        found_start = false
        for i in 1:chunks.num_rows
            y1 = chunk_tl[2] + (i-1) * chunks.chunk_height
            y2 = y1 + chunks.chunk_height
            if !found_start && y1 <= global_bbox.bottom_left[2] + global_bbox.height  <= y2
                found_start = true
                start_y_id = i
            end

            if found_start && y1 > global_bbox.bottom_left[2]
                end_y_id = i - 1
                break
            end
        end

        for col in start_x_id:end_x_id
            for row in start_y_id:end_y_id
                rect = RectangleCfg(
                    length=chunks.length,
                    height=chunks.height,
                    bottom_left=chunk_tl - col*up + (row-1)*right,
                )

                for (_, bbox, chunks_ids) in zip(bbox_vec, bbox_ids)
                    if check_intersection(rect, bbox)
                        push!(chunks_ids, CartesianIndex(row, col))
                    end
                end
            end
        end

        checker = ChunksChecker(chunks, chunks_ids)
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

function update_area_empty!(source::Source{T, C}) where {T, C<:PosChecker}
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

function update_area_empty!(source::Source{T, C}) where {T, C<:ChunksChecker}
    checker = source.checker 
    chunks = checker.chunks
    for (i, chunks_ids) in enumerate(checker.ids)
        is_empty = true
        source.is_empty[i] = is_empty
        for chunk_id in chunks_ids
            if chunks.num_particles_in_chunk[chunk_id] != 0
                is_empty = false
                break
            end
        end
        source.is_empty[i] = is_empty
    end
end

function process_source!(source::Nothing, state=nothing) end

function process_source!(source, state)
    update_area_empty!(source)
    for idx in eachindex(source.is_empty)
        if !source.is_empty[idx]
            continue
        end
        add_ring(state, source.bbox_spawn_pos[idx], get_spawn_pol(source.cfg.spawn_pol))
    end
end



end