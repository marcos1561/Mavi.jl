module MovableObjects

export LineState, MovableObject
export OverdampedCfg
export update_object!

using StaticArrays, LinearAlgebra
using Mavi.Configs

@kwdef struct OverdampedCfg{T}
    mu::T
end

abstract type Restriction end

abstract type MovableObjectState end

function update_object!(obj_state::MovableObjectState, dynamic_cfg::Nothing, restriction, force, int_cfg) end

mutable struct LineState{S, T} <: MovableObjectState
    p1::SVector{S, T}
    p2::SVector{S, T}
    normal::SVector{S, T}
    tangent::SVector{S, T}
    length::T
    vel::SVector{S, T}
end

function LineState(; p1, p2, vel) 
    D1, D2, D3 = length.((p1, p2, vel))
    if !(D1 == D2 == D3)
        error("p1, p2, and vel must have the same dimension")
    end
    D = D1

    T1, T2, T3 = eltype.((p1, p2, vel))
    T = promote_type(T1, T2, T3)

    p1 = SVector{D, T}(p1...)
    p2 = SVector{D, T}(p2...)
    vel = SVector{D, T}(vel...)
    
    PT = typeof(p1)
    line = LineState(p1, p2, zero(PT), zero(PT), zero(T), vel)

    update_state!(line)
    return line
end

function update_state!(line::LineState)
    dr = line.p2 - line.p1
    norm = sqrt(sum(abs2, dr))

    D = num_dimensions(line)
    if D == 2
        n = SVector(-dr[2], dr[1]) / norm
    elseif D == 3
        if dr[1] == 0 && dr[2] == 0
            other = SVector(0,1,0)
        else
            other = SVector(0,0,1)
        end
        n = cross(dr, other) / norm
    end
    
    t = dr / norm
    
    line.normal = n
    line.tangent = t
    line.length = norm
end

num_dimensions(line::LineState) = length(line.p1)
numerical_type(line::LineState) = eltype(line.p1)

function move_rigid!(line::LineState, dp)
    line.p1 += dp
    line.p2 += dp
end

function move!(line::LineState; dp1=nothing, dp2=nothing)
    if dp1 !== nothing
        line.p1 += dp1
    end
    if dp2 !== nothing
        line.p2 += dp2
    end
    update_state!(line)
end

function update_object!(line::LineState, dynamic_cfg::OverdampedCfg, restriction::Nothing, force, int_cfg)
    line.vel = dynamic_cfg.mu * force

    dt = int_cfg.dt
    cm = (line.p1 + line.p2) / 2

    new_cm = cm + dt * line.vel
    
    dr = new_cm - cm
    move_rigid!(line, dr)
    # line.p1 += dr
    # line.p2 += dr
end

function update_object!(line::LineState, dynamic_cfg::OverdampedCfg, restriction::Configs.Line2D, force, int_cfg)
    line.vel = dynamic_cfg.mu * force

    dt = int_cfg.dt
    cm = (line.p1 + line.p2) / 2

    rest_po = restriction.p1
    rest_tan = restriction.tangent

    new_cm = cm + dt * line.vel 
    new_cm = rest_po + dot((new_cm - rest_po), rest_tan) * rest_tan

    dr = new_cm - cm
    move_rigid!(line, dr)
    # line.p1 += dr
    # line.p2 += dr
end

function Configs.signed_pos(point, line::LineState)
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

mutable struct MovableObject{S, T, D, G<:MovableObjectState, P, R}
    force::SVector{S, T}
    external_force::SVector{S, T}
    dynamic_cfg::D
    state::G
    potential::P
    restriction::R
end
function MovableObject(; state, dynamic_cfg, potential, restriction, external_force=nothing)
    S, T = num_dimensions(state), numerical_type(state)
    force = SVector{S, T}(zeros(T, S))

    if external_force === nothing
        external_force = zero(force)
    else
        external_force = SVector{S, T}(external_force...)
    end

    MovableObject(force, external_force, dynamic_cfg, state, potential, restriction)
end

function update_object!(object::MovableObject, int_cfg)
    update_object!(object.state, object.dynamic_cfg, object.restriction, object.force, int_cfg)
end

end #MovableObjects