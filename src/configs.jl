"Collection of configurations"
module Configs

export SpaceCfg, DynamicCfg, IntegrationCfg

@kwdef mutable struct SpaceCfg
    length::Float64
    height::Float64
end

@kwdef struct DynamicCfg 
    ko::Float64
    ro::Float64
    ra::Float64
end

@kwdef struct IntegrationCfg 
    dt::Float64
end

end