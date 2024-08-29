"Collection of configurations"
module Configs

export SpaceCfg, DynamicCfg, AbstractIntCfg
export RectangleCfg, CircleCfg, HarmTruncCfg, LenJonesCfg
export IntCfg, ChunksIntCfg
export particle_radius, ChunksCfg

#
# Space Configurations 
#
abstract type SpaceCfg end

@kwdef struct RectangleCfg <: SpaceCfg
    length::Float64
    height::Float64
end

@kwdef struct CircleCfg <: SpaceCfg
    radius::Float64
end

"""Return bounding RectangleCfg for given space_cfg."""
function get_bounding_box(space_cfg::CircleCfg)
    r = space_cfg.radius * 2
    return RectangleCfg(r, r)
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

particle_radius(dynamic_cfg::HarmTruncCfg) = dynamic_cfg.ro/2
particle_radius(dynamic_cfg::LenJonesCfg) = dynamic_cfg.sigma/2

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

end