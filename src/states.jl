module States

using StaticArrays

export SecondLawState, SelfPropelledState, State

abstract type State{N, T} end

"Particles state for Newton's second law (positions and velocities)"
@kwdef struct SecondLawState{N, T} <: State{N, T}
    pos::Vector{SVector{N, T}}
    vel::Vector{SVector{N, T}}
end
function SecondLawState(pos::Matrix, vel::Matrix)
    @assert size(pos, 2) == size(vel, 2) "pos and vel must have the same number of particles"
    @assert size(pos, 1) == size(vel, 1) "pos and vel must have the same dimension"
    T = promote_type(eltype(pos), eltype(vel))

    pos = convert(Matrix{T}, pos)
    vel = convert(Matrix{T}, vel)

    SecondLawState(
        copy(reinterpret(SVector{size(pos, 1), T}, vec(pos))),
        copy(reinterpret(SVector{size(vel, 1), T}, vec(vel))),
    )
end

"Overdamped and self propelled state (positions and polarizations)"
@kwdef struct SelfPropelledState{N, T} <: State{N, T}
    pos::Vector{SVector{N, T}}
    pol_angle::Vector{T}
end
function SelfPropelledState(pos::Matrix, pos_angle)
    T = promote_type(eltype(pos), eltype(vel))

    pos = convert(Matrix{T}, pos)
    pos_angle = convert(Vector{T}, pos_angle)

    SelfPropelledState(
        copy(reinterpret(SVector{size(pos, 1), T}, vec(pos))),
        pos_angle,
    )
end

end