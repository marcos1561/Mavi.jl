"Collection of configurations"
module Configs

export GeometryCfg, DynamicCfg, AbstractIntCfg, GeometryCfg, WallType
export SpaceCfg, RectangleCfg, CircleCfg, LinesCfg 
export get_bounding_box, check_intersection, is_inside
export RigidWalls, PeriodicWalls, SlipperyWalls
export HarmTruncCfg, LenJonesCfg, SzaboCfg, RunTumbleCfg
export IntCfg, ChunksIntCfg, has_chunks
export particle_radius, ChunksCfg
export DeviceMode, Sequencial, Threaded

using StaticArrays

#
# Space Configurations 
#
abstract type GeometryCfg end

struct RectangleCfg{N, T<:Number} <: GeometryCfg
    length::T
    height::T
    bottom_left::SVector{N, T}
    size::SVector{N, T}
end
function RectangleCfg(;length, height, bottom_left=0, num_dims=2)
    if bottom_left === 0
        length, height = promote(length, height)
        bottom_left = SVector{num_dims, typeof(length)}(fill(0, num_dims))
    else
        length, height, bottom_left_correct_type... = promote(length, height, bottom_left...)
        num_dims = Base.length(bottom_left)
        bottom_left = SVector(bottom_left_correct_type)
    end

    RectangleCfg(length, height, bottom_left, SVector(length, height))
end

function Base.:+(a::RectangleCfg, b::RectangleCfg)
    max_x = max(a.bottom_left[1] + a.length, b.bottom_left[1] + b.length)
    max_y = max(a.bottom_left[2] + a.height, b.bottom_left[2] + b.height)
    min_x = min(a.bottom_left[1], b.bottom_left[1])
    min_y = min(a.bottom_left[2], b.bottom_left[2])

    height = max_y - min_y
    length = max_x - min_x
    RectangleCfg(length, height, SVector(min_x, min_y), SVector(length, height))
end

function check_intersection(r1::RectangleCfg, r2::RectangleCfg)
    intersect_x1 = r2.bottom_left[1] <= r1.bottom_left[1] <= r2.bottom_left[1] + r2.length
    intersect_x2 = r2.bottom_left[1] <= r1.bottom_left[1] + r1.length <= r2.bottom_left[1] + r2.length
    
    intersect_y1 = r2.bottom_left[2] <= r1.bottom_left[2] <= r2.bottom_left[2] + r2.height
    intersect_y2 = r2.bottom_left[2] <= r1.bottom_left[2] + r1.height <= r2.bottom_left[2] + r2.height

    return (intersect_x1 || intersect_x2) && (intersect_y1 || intersect_y2)
end

function is_inside(point, r::RectangleCfg; pad=0)
    is_x = r.bottom_left[1] - pad <= point[1] <= r.bottom_left[1] + r.length + pad
    is_y = r.bottom_left[2] - pad <= point[2] <= r.bottom_left[2] + r.height + pad
    return is_x && is_y
end

# struct Point2D{T}
#     x::T
#     y::T
# end
# Point2D(x, y) = Point2D(promote(x, y)...)
# Point2D(p) = Point2D(p...)

struct Line2D{T}
    p1::SVector{2, T}
    p2::SVector{2, T}
    normal::SVector{2, T}
    tangent::SVector{2, T}
    length::T
end
function Line2D(p1, p2) 
    if !(isa(p1, SVector))
        p1 = SVector(p1...)
    end
    if !(isa(p2, SVector))
        p2 = SVector(p2...)
    end
    
    # dx = p2.x - p1.x
    # dy = p2.y - p1.y
    # norm = sqrt(dx^2 + dy^2)
    dr = p2 - p1
    norm = sqrt(sum(abs2, dr))

    # n = [-dy, dx] / norm
    # t = [dx, dy] / norm

    n = SVector(-dr.y, dr.x) / norm
    t = dr / norm

    return Line2D(p1, p2, n, t, norm)
end

@kwdef struct LinesCfg{T} <: GeometryCfg
    lines::Vector{Line2D{T}}
    bbox::Union{RectangleCfg{2, T}, Nothing}
end
function LinesCfg(lines::AbstractVector; bbox=nothing)
    LinesCfg([Line2D(ps...) for ps in lines], bbox)
end

@kwdef struct CircleCfg <: GeometryCfg
    radius::Float64
end

"Return bounding RectangleCfg for given space_cfg."
function get_bounding_box(space_cfg::CircleCfg)
    r = space_cfg.radius
    length = 2 * r
    bottom_left = SVector(-r, -r)
    size_vec = SVector(length, length) 
    return RectangleCfg(length, length, bottom_left, size_vec)
end

function get_bounding_box(space_cfg::RectangleCfg{T}) where T
    return space_cfg
end

function get_bounding_box(space_cfg::LinesCfg)
    if !isnothing(space_cfg.bbox)
        return space_cfg.bbox
    end

    xs = []
    ys = []
    for line in space_cfg.lines
        push!(xs, line.p1.x, line.p2.x)
        push!(ys, line.p1.y, line.p2.y)
    end
    
    max_x, min_x = maximum(xs), minimum(xs)
    max_y, min_y = maximum(ys), minimum(ys)
    
    return RectangleCfg(
        length=max_x - min_x,
        height=max_y - min_y,
        bottom_left=(min_x, min_y),
    )
end

#
# Wall Type
#
abstract type WallType end

struct RigidWalls <: WallType end  
struct PeriodicWalls <: WallType end  
struct SlipperyWalls <: WallType end  

@kwdef struct SpaceCfg{W<:WallType, G<:GeometryCfg} 
    wall_type::W
    geometry_cfg::G
end

#
# Dynamic Configurations 
#
abstract type DynamicCfg end

"""
Truncated harmonic potencial configuration.

# Arguments
- ko:   
    Coupling constant

- ro:   
    Oscillation center

- ra:   
    Maximum distance for which the potential is nonzero
"""
@kwdef struct HarmTruncCfg{T<:Number} <: DynamicCfg 
    ko::T
    ro::T
    ra::T
end

"""
Lennard-Jones potential.

# Arguments
- sigma:   
    Distance between particles where the potential is zero.

- epsilon:   
    Depth of the potential well.
"""
@kwdef struct LenJonesCfg{T<:Number} <: DynamicCfg
    sigma::T
    epsilon::T
end

@kwdef struct SzaboCfg{T<:Number} <: DynamicCfg
    vo::T
    mobility::T
    relax_time::T
    k_rep::T
    k_adh::T
    r_eq::T
    r_max::T
    rot_diff::T
end

@kwdef struct RunTumbleCfg{T<:Number} <: DynamicCfg
    vo::T
    sigma::T
    epsilon::T
    tumble_rate::T
end

particle_radius(dynamic_cfg::HarmTruncCfg) = dynamic_cfg.ro/2
particle_radius(dynamic_cfg::LenJonesCfg) = dynamic_cfg.sigma * 2^(1/6) / 2
particle_radius(dynamic_cfg::SzaboCfg) = dynamic_cfg.r_eq/2
particle_radius(dynamic_cfg::RunTumbleCfg) = dynamic_cfg.sigma * 2^(1/6) / 2

#
# Integration Configuration
#
abstract type AbstractIntCfg end

abstract type DeviceMode end
struct Threaded <: DeviceMode end
struct Sequencial <: DeviceMode end

@kwdef struct ChunksCfg
    num_cols::Int
    num_rows::Int
end

@kwdef struct IntCfg{T<:Number, ChunkT<:Union{ChunksCfg, Nothing}, Device<:DeviceMode} <: AbstractIntCfg 
    dt::T
    chunks_cfg::ChunkT = nothing
    device::Device = Sequencial()
end

has_chunks(int_cfg::IntCfg{T, Nothing, D}) where {T, D} = false
has_chunks(int_cfg::IntCfg{T, ChunksCfg, D}) where {T, D} = true

end