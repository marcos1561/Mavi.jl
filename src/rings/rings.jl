module Rings
    
export RingsSystem, States, Configs, InitStates, Integration

import Mavi

@kwdef struct RingsInfo{T}
    continuos_pos::Array{T, 3}
    areas::Vector{Float64}
end

@kwdef struct DebugInfo{T}
    graph_pos::Matrix{T}
    graph_radius::Vector{Float64}
    graph_color::Vector{Symbol}
    graph_type::Vector{Int}
end

include("state.jl")
include("configs.jl")
include("integration.jl")
include("utils.jl")
include("init_state.jl")
include("view.jl")

function RingsSystem(;state, space_cfg, dynamic_cfg, int_cfg)
    info = RingsInfo(
        continuos_pos = similar(state.rings_pos),
        areas = Vector{Float64}(undef, size(state.rings_pos)[3]),
    )

    num_total_p = length(state.pos)
    debug_info = DebugInfo(
        graph_pos = similar(state.pos),
        graph_radius = Vector{Float64}(undef, num_total_p),
        graph_color = Vector{Symbol}(undef, num_total_p),
        graph_type = Vector{Int}(undef, num_total_p),
    )

    Mavi.System(
        state=state, 
        space_cfg=space_cfg, 
        dynamic_cfg=dynamic_cfg, 
        int_cfg=int_cfg, 
        info=info, 
        debug_info=debug_info,
    )
end

end