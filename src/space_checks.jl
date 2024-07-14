module SpaceChecks

using Mavi: State
using Mavi.Configs

export check_inside

"Return particles indices outside the given space configuration."
function outside_particles(state::State, space_cfg::RectangleCfg)
    x_out = @. (state.pos[1, :] < 0.0) | (state.pos[1, :] > space_cfg.length)
    y_out = @. (state.pos[2, :] < 0.0) | (state.pos[2, :] > space_cfg.height)
    out_ids = findall(x_out .| y_out)
    return out_ids
end

function outside_particles(state::State, space_cfg::CircleCfg)
    r2 = @. state.pos[1, :]^2 + state.pos[2, :]^2
    out_ids = findall(r2 .> space_cfg.radius^2)
    return out_ids
end

"""
Checks if all particles are inside the given space configuration. 

# Return
- `all_inside::bool`: True if all particles are inside, false otherwise.
- `out_ids::Vector{Int}`: Indices of outside particles.
"""
function check_inside(state::State, space_cfg::SpaceCfg)
    all_inside = true
    out_ids = outside_particles(state, space_cfg)
    if length(out_ids) != 0
        all_inside = false
    end
    return all_inside, out_ids
end

end