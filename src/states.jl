module States

export SecondLawState, SelfPropelledState, State

abstract type State{T} end

"Particles state for Newton's second law (positions and velocities)"
@kwdef struct SecondLawState{T} <: State{T}
    pos::Matrix{T}
    vel::Matrix{T}
end

"Overdamped and self propelled state (positions and polarizations)"
@kwdef struct SelfPropelledState{T} <: State{T}
    pos::Matrix{T}
    pol_angle::Vector{T}
end

end