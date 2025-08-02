module States

using StaticArrays

export SecondLawState, SelfPropelledState, State

abstract type State{N, T} end

"Particles state for Newton's second law (positions and velocities)"
@kwdef struct SecondLawState{N, T} <: State{N, T}
    pos::Vector{SVector{N, T}}
    vel::Vector{SVector{N, T}}
end

"Overdamped and self propelled state (positions and polarizations)"
@kwdef struct SelfPropelledState{N, T} <: State{N, T}
    pos::Matrix{T}
    pol_angle::Vector{T}
end

end