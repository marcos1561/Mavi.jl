"Collection of configurations"
module Configs

export GeometryCfg, DynamicCfg, AbstractIntCfg, GeometryCfg, WallType
export SpaceCfg, RectangleCfg, CircleCfg, LinesCfg 
export RigidWalls, PeriodicWalls, SlipperyWalls
export HarmTruncCfg, LenJonesCfg, SzaboCfg, RunTumbleCfg
export IntCfg, ChunksIntCfg, has_chunks
export particle_radius, ChunksCfg
export DeviceMode, Sequencial, Threaded

#
# Space Configurations 
#
abstract type GeometryCfg end

struct RectangleCfg{T<:Number} <: GeometryCfg
    length::T
    height::T
    bottom_left::Tuple{T, T}
end
function RectangleCfg(;length, height, bottom_left=(0, 0))
    args = promote(length, height, bottom_left...)
    RectangleCfg(args[1:end-2]..., args[end-1:end])
end


struct Point2D{T}
    x::T
    y::T
end
Point2D(x, y) = Point2D(promote(x, y)...)
Point2D(p) = Point2D(p...)

struct Line2D{T}
    p1::Point2D{T}
    p2::Point2D{T}
    normal::Point2D{T}
    tangent::Point2D{T}
    length::T
end
function Line2D(p1, p2) 
    if !(isa(p1, Point2D))
        p1 = Point2D(p1...)
    end
    if !(isa(p2, Point2D))
        p2 = Point2D(p2...)
    end
    
    dx = p2.x - p1.x
    dy = p2.y - p1.y
    norm = sqrt(dx^2 + dy^2)

    n = [-dy, dx] / norm
    t = [dx, dy] / norm
    return Line2D(p1, p2, Point2D(n), Point2D(t), norm)
end

@kwdef struct LinesCfg{T} <: GeometryCfg
    lines::Vector{Line2D{T}}
end
LinesCfg(points::AbstractVector) = LinesCfg([Line2D(ps...) for ps in points])


@kwdef struct CircleCfg <: GeometryCfg
    radius::Float64
end

"Return bounding RectangleCfg for given space_cfg."
function get_bounding_box(space_cfg::CircleCfg)
    r = space_cfg.radius
    length = 2 * r
    bottom_left = (-r, -r)
    return RectangleCfg(length, length, bottom_left)
end

function get_bounding_box(space_cfg::RectangleCfg{T}) where T
    return space_cfg
end

function get_bounding_box(space_cfg::LinesCfg)
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