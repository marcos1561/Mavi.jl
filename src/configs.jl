"Collection of configurations"
module Configs

export SpaceCfg, DynamicCfg, AbstracIntCfg
export RectangleCfg, CircleCfg, HarmTruncCfg, LenJonesCfg
export IntCfg, ChuncksIntCfg
export particle_radius, ChuncksCfg

#
# Space Configuraitons 
#
abstract type SpaceCfg end

@kwdef struct RectangleCfg <: SpaceCfg
    length::Float64
    height::Float64
end

@kwdef struct CircleCfg <: SpaceCfg
    radius::Float64
end

#
# Dynamic Configurations 
#
abstract type DynamicCfg end

"""
Truncated harmonic potencial configuration.

# Arguments
- ko: coupling constant
- ro: oscillation center
- ra: maximum distance for which the potential is nonzero
"""
@kwdef struct HarmTruncCfg <: DynamicCfg 
    ko::Float64
    ro::Float64
    ra::Float64
end

"""
Lennard-Jones potential

Arguments
- sigma: distance between particles where the potential is zero
- epsilon: depth of the potential well
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
abstract type AbstracIntCfg end

@kwdef struct IntCfg <: AbstracIntCfg 
    dt::Float64
end

@kwdef struct ChuncksCfg
    num_cols::Int
    num_rows::Int
end

@kwdef struct ChuncksIntCfg <: AbstracIntCfg
    dt::Float64
    chuncks_cfg::ChuncksCfg
end

end