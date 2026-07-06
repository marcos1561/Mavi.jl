"Collection of configurations"
module Configs

export GeometryCfg, DynamicCfg, AbstractIntCfg, GeometryCfg, WallType
export SpaceCfg, RectangleCfg, CircleCfg, LinesCfg, ManyGeometries 
export get_bounding_box, check_intersection, signed_pos, is_inside
export RigidWalls, PeriodicWalls, SlipperyWalls, ManyWalls
export ForceWalls, PotentialWalls, process_dist
export get_main_wall, get_main_geometry
export get_space_data, SpaceData
export get_potential_cfg

export HarmTruncCfg, LenJonesCfg, SzaboCfg, RunTumbleCfg
export IntCfg, ChunksCfg, has_chunks, get_chunks_cfg_from_cell_size
export DeviceMode, Sequencial, Threaded
export particle_radius, potential_force, maximum_interaction_distance

using StaticArrays, StructTypes, Random
using Mavi.States

# =
# Geometries 
# =
abstract type SpaceData end

get_space_data(space_cfg) = nothing  
    
abstract type GeometryCfg end

function is_inside(point, geometry::GeometryCfg; pad) 
    throw(NotImplementedError("is_inside not implemented for $(typeof(geometry))"))
end

function num_dimensions(geometry::GeometryCfg) 
    throw(NotImplementedError("num_dimensions not implemented for $(typeof(geometry))"))
end

function check_intersection(geometry_1::GeometryCfg, geometry_2::GeometryCfg)
    throw(NotImplementedError("check_intersection not implemented for $(typeof(geometry_1))-$(typeof(geometry_2))"))
end

function get_bounding_box(geometry::GeometryCfg)
    throw(NotImplementedError("get_bounding_box not implemented for $(typeof(geometry))"))
end

function signed_pos(point, geometry::GeometryCfg)
    throw(NotImplementedError("signed_pos not implemented for $(typeof(geometry_1))"))
end

function random_point(geometry::GeometryCfg)
    throw(NotImplementedError("random_point not implemented for $(typeof(geometry))"))
end

function get_numerical_type(geometry::GeometryCfg) 
    throw(NotImplementedError("get_numerical_type not implemented for $(typeof(geometry))"))
end 

function get_num_dimensions(geometry::GeometryCfg) 
    throw(NotImplementedError("get_num_dimensions not implemented for $(typeof(geometry))"))
end
    
function random_points(geometry::GeometryCfg, n)
    T, D = get_numerical_type(geometry), get_num_dimensions(geometry)
    points = Vector{SVector{D, T}}(undef, n)
    for i in 1:n
        points[i] = random_point(geometry)
    end
    return points
end

function random_points_no_overlap(geometry::GeometryCfg, n; radius=0, max_attempts=10, offset=0)
    T, D = get_numerical_type(geometry), get_num_dimensions(geometry)
    points = Vector{SVector{D, T}}(undef, n)
    r2 = radius^2
    for i in 1:n
        found_point = false
        num_attempts = 0
        while !found_point
            if num_attempts >= max_attempts
                error("Could not find a valid point!")
            end
            num_attempts += 1

            found_point = true
            p_try = random_point(geometry)
            for p in points[1:(i-1)]
                if sum((p_try - p).^2) < r2
                    found_point = false
                    break
                end
            end
            
            if signed_pos(p_try, geometry)[2] < offset
                found_point = false
            else
                points[i] = p_try
            end
        end
    end
    return points
end

struct ManyGeometries{G <: Tuple} <: GeometryCfg
    list::G
end

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

function RectangleCfg(points; offset=0)
    x = [points[i][1] for i in eachindex(points)]
    y = [points[i][2] for i in eachindex(points)]

    max_x = maximum(x) + offset
    min_x = minimum(x) - offset
    max_y = maximum(y) + offset
    min_y = minimum(y) - offset
    RectangleCfg(
        length=max_x-min_x, 
        height=max_y-min_y, 
        bottom_left=(min_x, min_y),
    )
end

get_numerical_type(geometry::RectangleCfg{N, T}) where {N, T} = T 
get_num_dimensions(geometry::RectangleCfg{N, T}) where {N, T} = N 

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

function random_point(geometry_cfg::RectangleCfg{N, T}) where {N, T}
    geometry_cfg.bottom_left + rand(SVector{N, T}) .* geometry_cfg.size
end

function signed_pos(point, geometry_cfg::RectangleCfg)
    dr = point - geometry_cfg.bottom_left

    dist_faces_1 = -dr 
    dist_faces_2 = dr - geometry_cfg.size
    dist_faces = vcat(dist_faces_1, dist_faces_2)

    min_dist, idx = findmin(abs.(dist_faces))
    sign_dist = dist_faces[idx] 

    D = length(dr)
    dr = zeros(eltype(dr), D)
    face_type = idx <= D ? 1 : 2
    possible_signs = [-1, 1]
    dr[mod1(idx, D)] = possible_signs[face_type] * sign_dist

    return dr, min_dist, sign(sign_dist)
end

abstract type AbstractLine end

struct Line2D{T} <: AbstractLine
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
    
    dr = p2 - p1
    norm = sqrt(sum(abs2, dr))

    n = SVector(-dr.y, dr.x) / norm
    t = dr / norm

    return Line2D(p1, p2, n, t, norm)
end

function signed_pos(point, line::Line2D)
    p1 = line.p1
            
    dr = point - p1
    delta_t = sum(dr .* line.tangent)
    
    if delta_t > 0
        if delta_t < line.length
            base_pos = p1 + line.tangent * delta_t
            dr = point - base_pos
        else
            dr = point - line.p2
        end
    end

    return dr, sqrt(sum(dr.^2)), 1
end

@kwdef struct LinesCfg{T} <: GeometryCfg
    lines::Vector{Line2D{T}}
    bbox::Union{RectangleCfg{2, T}, Nothing}
end
function LinesCfg(lines::AbstractVector; bbox=nothing)
    LinesCfg([Line2D(ps...) for ps in lines], bbox)
end

struct CircleCfg{N, T<:Number} <: GeometryCfg
    radius::T
    center::SVector{N, T}
end
function CircleCfg(;radius, center)
    radius, center... = promote(radius, center...)
    CircleCfg(radius, SVector(center...))
end

function signed_pos(point, geometry_cfg::CircleCfg)
    dr = point - geometry_cfg.center
    dist = sum(dr.^2)^.5
    dr_hat = dr / dist

    dr = dr - dr_hat * geometry_cfg.radius
    sign_dist = dist - geometry_cfg.radius
    
    return dr, abs(sign_dist), sign(sign_dist)
end

function is_inside(point, r::CircleCfg; pad=0)
    dist_sqr = sum((point - r.center).^2)    
    return dist_sqr <= (r.radius + pad)^2 
end

check_intersection(r1::RectangleCfg, r2::CircleCfg) = check_intersection(r2, r1)
function check_intersection(r1::CircleCfg, r2::RectangleCfg)
    # Clamp the circle's center to the rectangle
    closest_p = clamp.(r1.center, r2.bottom_left, r2.bottom_left + r2.size)

    # Calculate the distance between the circle's center and the closest point
    distance_sqr = sum((closest_p - r1.center).^2) 

    # Check if the distance is less than or equal to the circle's radius
    return distance_sqr <= r1.radius^2
end

function check_intersection(r1::CircleCfg, r2::CircleCfg)
    diff = r1.center - r2.center
    dist_sqr = sum(diff.^2)
    return dist_sqr < (r1.radius + r2.radius)^2
end

"Return bounding RectangleCfg for given space_cfg."
function get_bounding_box(geometry_cfg::CircleCfg)
    r = geometry_cfg.radius
    length = 2 * r
    bottom_left = geometry_cfg.center - SVector(r, r)
    size_vec = SVector(length, length) 
    return RectangleCfg(length, length, bottom_left, size_vec)
end

function get_bounding_box(geometry_cfg::RectangleCfg{T}) where T
    return geometry_cfg
end

function get_bounding_box(geometry_cfg::LinesCfg)
    if !isnothing(geometry_cfg.bbox)
        return geometry_cfg.bbox
    end

    xs = []
    ys = []
    for line in geometry_cfg.lines
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

function get_bounding_box(geometry_cfg::ManyGeometries)
    bbox = get_bounding_box(geometry_cfg.list[1])
    for g in geometry_cfg.list[2:end]
        bbox += get_bounding_box(g)
    end
    return bbox
end

# =
# Wall Type
# =

abstract type WallType end
StructTypes.StructType(::Type{W}) where W <: WallType = StructTypes.Struct()

struct RigidWalls <: WallType end  
struct PeriodicWalls <: WallType end  
struct SlipperyWalls <: WallType end  

struct ManyWalls{W <: Tuple} <: WallType 
    list::W
end  

# =
# Force Walls
# =

abstract type ForceWalls <: WallType end

abstract type PotentialWallMode end
StructTypes.StructType(::Type{P}) where P <: PotentialWallMode = StructTypes.Struct()

struct Outside <: PotentialWallMode end
struct Inside <: PotentialWallMode end
struct Repulsion <: PotentialWallMode end

process_dist(mode::Outside, dist, inside_flag) = inside_flag * dist
process_dist(mode::Inside, dist, inside_flag) = -inside_flag * dist
process_dist(mode::Repulsion, dist, inside_flag) = dist

struct PotentialWalls{P, M<:PotentialWallMode} <: ForceWalls 
    potential::P
    mode::M
end
function PotentialWalls(;potential, mode=:repulsion)
    if mode isa Symbol
        mode_map = (
            outside=Outside(),
            inside=Inside(),
            repulsion=Repulsion(),
        )
        
        mode = mode_map[mode]
    end

    PotentialWalls(potential, mode)
end
    
# =
# Space Configuration
# = 

@kwdef struct SpaceCfg{W<:WallType, G<:GeometryCfg} 
    wall_type::W
    geometry_cfg::G
end
function SpaceCfg(spaces_cfg_list)
    walls = []
    geometries = []

    for (w, g) in spaces_cfg_list
        push!(walls, w)
        push!(geometries, g)
    end

    SpaceCfg(
        wall_type=ManyWalls(Tuple(walls)),
        geometry_cfg=ManyGeometries(Tuple(geometries)),
    )
end

"Wall type used to calculate things as chunks neighbors."
get_main_wall(wall::WallType) = wall
get_main_wall(wall::ManyWalls) = wall.list[1]
get_main_wall(space_cfg::SpaceCfg) = get_main_wall(space_cfg.wall_type)

"Geometry used as border."
get_main_geometry(geometry_cfg::GeometryCfg) = geometry_cfg
get_main_geometry(geometry_cfg::ManyGeometries) = geometry_cfg.list[1]
get_main_geometry(space_cfg::SpaceCfg) = get_main_geometry(space_cfg.geometry_cfg)

# ==
# Dynamic Configurations 
# =

abstract type DynamicCfg end
maximum_interaction_distance(potential::DynamicCfg) = @error "Not Implemented for $(typeof(potential))"

abstract type PotentialCfg <: DynamicCfg end
potential_force(dr, dist, potential::PotentialCfg, system, p_id) = potential_force(dr, dist, potential)

function potential_force(dr, potential::PotentialCfg)
    potential_force(dr, sqrt(sum(dr.^2)), potential)
end

"""
Truncated harmonic potential configuration.

# Arguments
- k_rep:   
    Repulsive force constant

- k_atr:   
    Attraction force constant

- dist_eq:   
    Equilibrium distance (where the force is zero)

- dist_max:   
    Distance beyond which the potential is truncated
"""
struct HarmTruncCfg{T<:Number} <: PotentialCfg
    k_rep::T
    k_atr::T
    dist_eq::T
    dist_max::T
end
function HarmTruncCfg(;k_rep, k_atr, dist_eq, dist_max)
    HarmTruncCfg(promote(k_rep, k_atr, dist_eq, dist_max)...)
end
function HarmTruncCfg{T}(; k_rep, k_atr, dist_eq, dist_max) where {T<:Number}
    HarmTruncCfg{T}(T(k_rep), T(k_atr), T(dist_eq), T(dist_max))
end

function potential_force(dr, dist, potential::HarmTruncCfg)
    if dist > potential.dist_max
        return zero(dr)
    end

    dist_eq = potential.dist_eq

    if dist < dist_eq
        fmod = -potential.k_rep * (dist/dist_eq - 1)
    else
        fmod = -potential.k_atr * (dist/dist_eq - 1)
    end

    return fmod / dist * dr
end

"""
Lennard-Jones potential.

# Arguments
- sigma:   
    Distance between particles where the potential is zero.

- epsilon:   
    Depth of the potential well.
"""
struct LenJonesCfg{T<:Number} <: PotentialCfg
    sigma::T
    epsilon::T
    sigma6::T
    sigma12::T
end
function LenJonesCfg(sigma, epsilon)
    sigma, epsilon = promote(sigma, epsilon)
    LenJonesCfg(sigma, epsilon, sigma^6, sigma^12)
end
LenJonesCfg(;sigma, epsilon) = LenJonesCfg(sigma, epsilon)

function potential_force(dr, dist, potential::LenJonesCfg)
    sigma6 = potential.sigma6
    sigma12 = potential.sigma12
    epsilon = potential.epsilon

    d2 = dist * dist
    d4 = d2 * d2
    d6 = d4 * d2
    d7 = d6 * dist
    d13 = d6 * d7

    # Force modulus
    fmod = 4*epsilon*(12*sigma12/d13 - 6*sigma6/d7)

    return fmod/dist * dr 
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

maximum_interaction_distance(p::HarmTruncCfg) = p.dist_max
maximum_interaction_distance(p::LenJonesCfg) = Inf
maximum_interaction_distance(p::SzaboCfg) = p.r_max

particle_radius(dynamic_cfg::HarmTruncCfg) = dynamic_cfg.dist_eq / 2.0
particle_radius(dynamic_cfg::LenJonesCfg) = dynamic_cfg.sigma * 2^(1/6) / 2
particle_radius(dynamic_cfg::SzaboCfg) = dynamic_cfg.r_eq/2
particle_radius(dynamic_cfg::RunTumbleCfg) = dynamic_cfg.sigma * 2^(1/6) / 2

# ==
# Potential Finder
# ==
abstract type PotentialFinder{T} <: DynamicCfg end

struct PotentialMatrix{T} <: PotentialFinder{T} 
    matrix::Matrix{T}
end

list_potentials(potentials::PotentialCfg) = [potentials]

function list_potentials(potentials::PotentialMatrix)
    mat = potentials.matrix
    return [mat[i, j] for i in 1:size(mat, 1) for j in i:size(mat, 2)]
end

list_self_potentials(potentials::PotentialCfg) = [potentials]

function list_self_potentials(potentials::PotentialMatrix)
    mat = potentials.matrix
    return [mat[i, i] for i in axes(mat, 1)]
end

@inline function get_potential_cfg(type_1, type_2, potential_finder::PotentialMatrix)
    potential_finder.matrix[type_1, type_2]
end

@inline get_potential_cfg(type_1, type_2, potential::PotentialCfg) = potential
@inline get_potential_cfg(potential, state, pid) = potential


struct PotentialVector{T} <: PotentialFinder{T} 
    vector::Vector{T}
end

list_potentials(potentials::PotentialVector) = potentials.vector

function get_potential_cfg(potentials::PotentialVector, state, pid) 
    p_type = get_particle_type(state, pid)
    return potentials.vector[p_type]
end 


# =
# Integration Configuration
# =
abstract type AbstractIntCfg end

abstract type DeviceMode end
StructTypes.StructType(::Type{D}) where D <: DeviceMode = StructTypes.Struct()

struct Threaded <: DeviceMode end
struct Sequencial <: DeviceMode end

@kwdef struct ChunksCfg
    num_cols::Int
    num_rows::Int
end
function get_chunks_cfg_from_cell_size(cell_size, rect_cfg::RectangleCfg)
    num_cells = floor.(Int, rect_cfg.size ./ cell_size)
    return ChunksCfg(num_cells...)
end

@kwdef struct IntCfg{T<:Number, ChunkT<:Union{ChunksCfg, Nothing}, Device<:DeviceMode, ExtraT} <: AbstractIntCfg 
    dt::T
    chunks_cfg::ChunkT = nothing
    device::Device = Sequencial()
    extra::ExtraT = nothing
end

has_chunks(int_cfg::IntCfg{T, Nothing, D, E}) where {T, D, E} = false
has_chunks(int_cfg::IntCfg{T, ChunksCfg, D, E}) where {T, D, E} = true

end # Configs