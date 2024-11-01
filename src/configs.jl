"Collection of configurations"
module Configs

export GeometryCfg, DynamicCfg, AbstractIntCfg, GeometryCfg, WallType
export SpaceCfg, RectangleCfg, CircleCfg, RigidWalls, PeriodicWalls
export HarmTruncCfg, LenJonesCfg, SzaboCfg, RunTumbleCfg
export IntCfg, ChunksIntCfg, has_chunks
export particle_radius, ChunksCfg

#
# Space Configurations 
#
abstract type GeometryCfg end

@kwdef struct RectangleCfg <: GeometryCfg
    length::Float64
    height::Float64
    bottom_left::Tuple{Float64, Float64}=(0.0, 0.0)
end

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

function get_bounding_box(space_cfg::RectangleCfg)
    return space_cfg
end

#
# Wall Type
#
abstract type WallType end

struct RigidWalls <: WallType end  
struct PeriodicWalls <: WallType end  

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
@kwdef struct HarmTruncCfg <: DynamicCfg 
    ko::Float64
    ro::Float64
    ra::Float64
end

"""
Lennard-Jones potential.

# Arguments
- sigma:   
    Distance between particles where the potential is zero.

- epsilon:   
    Depth of the potential well.
"""
@kwdef struct LenJonesCfg <: DynamicCfg
    sigma::Float64
    epsilon::Float64
end

@kwdef struct SzaboCfg <: DynamicCfg
    vo::Float64
    mobility::Float64
    relax_time::Float64
    k_rep::Float64
    k_adh::Float64
    r_eq::Float64
    r_max::Float64
    rot_diff::Float64
end

@kwdef struct RunTumbleCfg <: DynamicCfg
    vo::Float64
    sigma::Float64
    epsilon::Float64
    tumble_rate::Float64
end

particle_radius(dynamic_cfg::HarmTruncCfg) = dynamic_cfg.ro/2
particle_radius(dynamic_cfg::LenJonesCfg) = dynamic_cfg.sigma * 2^(1/6) / 2
particle_radius(dynamic_cfg::SzaboCfg) = dynamic_cfg.r_eq/2
particle_radius(dynamic_cfg::RunTumbleCfg) = dynamic_cfg.sigma * 2^(1/6) / 2

#
# Integration Configuration
#
abstract type AbstractIntCfg end

@kwdef struct IntCfg <: AbstractIntCfg 
    dt::Float64
end

@kwdef struct ChunksCfg
    num_cols::Int
    num_rows::Int
end

@kwdef struct ChunksIntCfg <: AbstractIntCfg
    dt::Float64
    chunks_cfg::ChunksCfg
end

has_chunks(int_cfg::IntCfg) = false
has_chunks(int_cfg::ChunksIntCfg) = true

end