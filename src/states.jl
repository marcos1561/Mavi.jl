module States

export StaticArrays

using StaticArrays

export SecondLawState, SelfPropelledState, State
export ActiveState, get_particles_ids, get_num_total_particles, update_ids!

# = 
# Active State used by users
# =

struct ActiveState{M}
    mask::M
end
ActiveState() = ActiveState(true)

get_active_mask(active_state::ActiveState{Bool}, num_entities) = fill(active_state.mask, num_entities)
get_active_mask(active_state::ActiveState{T}, num_entities=nothing) where T <: AbstractVector = vec(active_state.mask)

# = 
# Particles Ids for variable number of particles
# = 

abstract type AbstractParticleIds end

@kwdef mutable struct ParticleIds <: AbstractParticleIds
    mask::Vector{Bool}
    ids::Vector{Int}
    num::Int
end

get_ids(part_ids) = part_ids
function get_ids(part_ids::ParticleIds)
    ids = @view part_ids.ids[1:part_ids.num]
    return ids
end

get_num(part_ids) = length(part_ids)
get_num(part_ids::ParticleIds) = part_ids.num

function update_part_ids!(part_ids) end 
function update_part_ids!(part_ids::ParticleIds) 
    idx = 1
    for (p_id, is_active) in enumerate(part_ids.mask)
        if is_active
            part_ids.ids[idx] = p_id
            idx += 1
        end 
    end
    part_ids.num = idx - 1
end

function get_particles_ids_obj(active_state, pos)
    if isnothing(active_state)
        part_ids = eachindex(pos)
    else
        num_p = length(pos)

        part_ids = ParticleIds(
            mask=get_active_mask(active_state, num_p),
            ids=Vector{Int}(undef, num_p),
            num=0
        )
        update_part_ids!(part_ids)
    end

    return part_ids
end

abstract type State{N, T} end

"Particles state for Newton's second law (positions and velocities)"
struct SecondLawState{N, T, PID <: Union{ParticleIds, AbstractVector}} <: State{N, T}
    pos::Vector{SVector{N, T}}
    vel::Vector{SVector{N, T}}
    part_ids::PID
end
function SecondLawState(pos, vel, active_state=nothing)
    SecondLawState(pos, vel, get_particles_ids_obj(active_state, pos))
end
function SecondLawState(pos::Matrix, vel::Matrix, active_state=nothing)
    @assert size(pos, 2) == size(vel, 2) "pos and vel must have the same number of particles"
    @assert size(pos, 1) == size(vel, 1) "pos and vel must have the same dimension"
    
    T = promote_type(eltype(pos), eltype(vel))

    if T <: Int
        T = Float64
    end

    pos = convert(Matrix{T}, pos)
    vel = convert(Matrix{T}, vel)

    pos = copy(reinterpret(SVector{size(pos, 1), T}, vec(pos)))
    vel = copy(reinterpret(SVector{size(vel, 1), T}, vec(vel)))   

    SecondLawState(pos, vel, active_state)
end
SecondLawState(; pos, vel, active_state=nothing) = SecondLawState(pos, vel, active_state)

"Overdamped and self propelled state (positions and polarizations)"
struct SelfPropelledState{N, T, PID <: Union{ParticleIds, AbstractVector}} <: State{N, T}
    pos::Vector{SVector{N, T}}
    pol_angle::Vector{T}
    part_ids::PID
end
function SelfPropelledState(;pos, pol_angle, active_state=nothing)
    if pos isa Matrix
        pos = copy(reinterpret(SVector{size(pos, 1), eltype(pos)}, vec(pos)))
    end

    T = promote_type(eltype(pos[1]), eltype(pol_angle))
    N = length(pos[1])
    
    pos = map(p -> convert(SVector{N, T}, p), pos)
    pol_angle = convert(Vector{T}, pol_angle)

    SelfPropelledState(
        pos,
        pol_angle,
        get_particles_ids_obj(active_state, pos),
    )
end

update_ids!(state) = update_part_ids!(state.part_ids)

get_particles_ids(state) = get_ids(state.part_ids)
get_num_total_particles(state) = get_num(state.part_ids)

end