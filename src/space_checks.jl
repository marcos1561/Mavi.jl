module SpaceChecks

using Mavi.States: State
using Mavi.Configs

export check_inside

"Return particles indices outside the given geometry configuration."
function outside_particles(state::State, geometry_cfg::RectangleCfg) 
    bottom_left = geometry_cfg.bottom_left
    top_right = bottom_left .+ (geometry_cfg.length, geometry_cfg.height)
    out_ids::Vector{Int} = []
    for idx in eachindex(state.pos)
        pos = state.pos[idx]
        if any(pos .< bottom_left) || any(pos .> top_right)
            push!(out_ids, idx)
        end
    end

    # x_out = @. (state.pos[1, :] - geometry_cfg.bottom_left[1] < 0.0) | (state.pos[1, :] - geometry_cfg.bottom_left[1] > geometry_cfg.length)
    # y_out = @. (state.pos[2, :] - geometry_cfg.bottom_left[2] < 0.0) | (state.pos[2, :] - geometry_cfg.bottom_left[2] > geometry_cfg.height)
    # out_ids = findall(x_out .| y_out)
    
    return out_ids
end

function outside_particles(state::State, geometry_cfg::CircleCfg)
    r2 = @. state.pos[1, :]^2 + state.pos[2, :]^2
    out_ids = findall(r2 .> geometry_cfg.radius^2)
    return out_ids
end

function outside_particles(state::State, geometry_cfg)
    return []
end

"""
Checks if all particles are inside the given geometry configuration. 

# Return
- all_inside::bool
    True if all particles are inside, false otherwise.

- out_ids::Vector{Int} 
    Indices of outside particles.
"""
function check_inside(state::State, geometry_cfg::GeometryCfg)
    all_inside = true
    out_ids = outside_particles(state, geometry_cfg)
    if length(out_ids) != 0
        all_inside = false
    end
    return all_inside, out_ids
end

end