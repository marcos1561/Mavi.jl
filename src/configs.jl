"Collection of configurations"
module Configs

export SpaceCfg, DynamicCfg, IntegrationCfg
export RectangleCfg, CircleCfg, HarmTruncCfg
export particle_radius

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
# Dynamic Configuraitons 
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
particle_radius(dynamic_cfg::HarmTruncCfg) = dynamic_cfg.ro/2 

#
# Integration Configuration
#
@kwdef struct IntegrationCfg 
    dt::Float64
end

end