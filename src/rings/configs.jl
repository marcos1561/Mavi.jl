module Configs

export RingsIntCfg, RingsCfg, has_types_cfg, get_spring_pars
export get_area0, get_equilibrium_p0, get_particle_radius
export get_equilibrium_area, get_particles_area_contribution, get_ring_radius
export HarmTruncCfg, PotentialMatrix
export IntCfg, InvasionsCfg

using NLsolve, StructTypes, JSON3

import Mavi.Configs as MaviCfg

import Mavi.Rings.States: ring_num_particles

using Reexport
@reexport using Mavi.Configs

# =
# Rings Configs
# = 

struct RingsCfg{U<:Union{AbstractVector, Number}, T<:MaviCfg.PotentialCfg, PF<:Union{MaviCfg.PotentialFinder{T}, T}} <: DynamicCfg
    p0::U
    relax_time::U
    vo::U
    mobility::U
    rot_diff::U
    k_area::U
    k_spring::U
    l_spring::U
    num_particles::Union{Int, Vector{Int}}
    num_types::Int
    interaction_finder::PF
end

"""
Rings configurations.

# Arguments
- p0:
    Parameter controlling how flexible the ring membrane is. It is define as
    
    p0 = n * l / A^(1/2) 
    
    where n is the number of particles, l the spring length and A the ring
    area (viewing the ring as a polygon with vertices at its particles positions).
    
    It has an equilibrium value p_eq (see docs for `get_equilibrium_p0`), if
    - p0 > p_eq: The ring is inflated and its membrane is not very flexible.
    - p0 < p_eq: The ring can be deformed without resistance, therefore its membrane is very flexible.

- relax_time:
    Relaxation time between the polarization and ring velocity. Lower values
    induce more collective motion.  

- vo:
    Activity parameter: specifically, `vo` is the magnitude of the ring's velocity.

- mobility:
    Factor controlling how important the forces are to dynamics.

- rot_diff:
    Noise intensity applied to the polarization angle (θ). The noise
    is added to the time derivative of θ, as a term of the form 

    sqrt(2 * rot_diff) * ξ

    where ξ is a gaussian white noise.

- k_area: 
    Ring area stiffness. This parameters enters in the potential which
    tends to preserve the ring area
    
    k_area/2 * (A - A_0)^2

    where A is the ring area and A_0 its equilibrium value.

- k_spring:
    Stiffness of the springs that connect adjacent particles in a ring.

- l_spring:
    Equilibrium length of the springs that connect adjacent particles in a ring.
    
- interaction_finder:
    An object defining the interaction between different rings. It can be
    
    - A potential configuration: It this case, this potential will be used by every
    interaction between rings.
    
    - A potential finder: This object is used when there are different types of rings
    that interact differently depending on the type. This object defines the potential
    configuration given the types interacting. Currently, the only potential finder available
    is the `InteractionMatrix`, which is a matrix of potential configuration, where the element
    with index (i, j) has the potential configuration for an interaction between types i and j.

- num_particles:
    How many particles a ring has. If there are more than one type of ring, then a vector
    must be provided, where the i-th element is the number of particles for the i-th ring type.
    If not given, this parameter will be automatically inferred by the number of particles in `state`.
        
- NumType:
    Numerical type used. By default, it is automatically calculated using the parameters types. 
"""
function RingsCfg(;
    p0, relax_time, vo, mobility, rot_diff, k_area, k_spring, l_spring, 
    interaction_finder, num_particles=-1, NumType=nothing,
    )
    args = (p0, relax_time, vo, mobility, rot_diff, k_area, k_spring, l_spring)
    
    num_types = maximum([length(arg) for arg in args])
    if num_types > 1 && !all(x -> x isa AbstractVector, args)
        new_args = []
        for arg in args
            if arg isa Number
                arg = [arg for _ in 1:num_types]
            end
            push!(new_args, arg)
        end
        args = new_args
    end

    U = typeof(args[1])
    
    if num_particles != -1
        if U <: Number && !(num_particles isa Int)
            throw(ArgumentError("If U is Number, num_particles must be Int"))
        elseif U <: AbstractVector && !(num_particles isa AbstractVector)
            throw(ArgumentError("If U is AbstractVector, num_particles must be Vector"))
        elseif num_types != length(num_particles)
            np = length(num_particles)
            throw(ArgumentError("`length(num_particles)=$np`, but there exists $num_types types"))
        end
    end

    if isnothing(NumType)
        NumType = promote_type(map(eltype, args)...)
    end
    
    if U <: AbstractVector
        args = [convert(Vector{NumType}, arg) for arg in args]
    else
        args = [convert(NumType, arg) for arg in args]
    end

    RingsCfg(args..., num_particles, num_types, interaction_finder)
end
function RingsCfg(cfg::RingsCfg, num_particles)
    RingsCfg(
        p0=cfg.p0,
        relax_time=cfg.relax_time,
        vo=cfg.vo,
        mobility=cfg.mobility,
        rot_diff=cfg.rot_diff,
        k_area=cfg.k_area,
        k_spring=cfg.k_spring,
        l_spring=cfg.l_spring,
        interaction_finder=cfg.interaction_finder,
        num_particles=num_particles,
    )
end

function MaviCfg.maximum_interaction_distance(dynamic_cfg::RingsCfg)
    return maximum(
        p -> maximum_interaction_distance(p),
        list_potentials(dynamic_cfg.interaction_finder)
    )
end

@inline has_types_cfg(dynamic_cfg::RingsCfg{U, T, F}) where {U<:Number, T, F} = false
@inline has_types_cfg(dynamic_cfg::RingsCfg{U, T, F}) where {U<:AbstractVector, T, F} = true

@inline function get_ring_prop_by_name(dynamic_cfg::RingsCfg{U, T, I}, name, type=nothing) where {U<:Number, T, I}
    return getfield(dynamic_cfg, name)
end
@inline function get_ring_prop_by_name(dynamic_cfg::RingsCfg{U, T, I}, name, type) where {U<:AbstractVector, T, I}
    return getfield(dynamic_cfg, name)[type]
end

function get_rings_property(dynamic_cfg::RingsCfg{U, T, I}, prop_func; dtype=nothing) where {U<:Number, T, I}
    if isnothing(dtype)
        dtype = U
    end
    dtype(prop_func(dynamic_cfg, nothing))
end

function get_rings_property(dynamic_cfg::RingsCfg{U, T, I}, prop_func; dtype=nothing) where {U, T, I}
    if isnothing(dtype)
        if U <: Number
            dtype = U
        else 
            dtype = eltype(U)
        end
    end

    prop_vec = Vector{dtype}(undef, dynamic_cfg.num_types)
    for i in 1:dynamic_cfg.num_types
        prop_vec[i] = prop_func(dynamic_cfg, i)
    end
    
    if length(prop_vec) == 1
        prop_vec = prop_vec[1]
    end

    return prop_vec
end

function get_spring_pars(dynamic_cfg::RingsCfg{U, T, F}, state=nothing, ring_id=nothing) where {U<:Number, T, F}
    return dynamic_cfg.k_spring, dynamic_cfg.l_spring
end

function get_spring_pars(dynamic_cfg::RingsCfg{U, T, F}, state, ring_id) where {U<:AbstractVector, T, F}
    t = state.types[ring_id]
    return dynamic_cfg.k_spring[t], dynamic_cfg.l_spring[t]
end


@inline function get_particle_radius(dynamic_cfg::RingsCfg, type=nothing)
    inter = get_potential_cfg(type, type, dynamic_cfg.interaction_finder)
    return particle_radius(inter)
end

MaviCfg.particle_radius(dynamic_cfg::RingsCfg) = get_rings_property(dynamic_cfg, get_particle_radius)
MaviCfg.entity_radius(dynamic_cfg::RingsCfg) = get_rings_property(dynamic_cfg, get_ring_radius)


ring_num_particles(dynamic_cfg::RingsCfg, type) = ring_num_particles(dynamic_cfg.num_particles, type)
ring_num_particles(dynamic_cfg::RingsCfg) = get_rings_property(dynamic_cfg, ring_num_particles, dtype=Int)

"""
Returns the equilibrium area of the area force for
the given number of particles, taking into account `p0`.
"""
get_area0(num_particles, l_spring, p0) = (num_particles * l_spring / p0)^2

function get_area0(dynamic_cfg::RingsCfg{U, T, I}, type=nothing) where {U, T, I}
    num_particles = ring_num_particles(dynamic_cfg, type)
    l_spring = get_ring_prop_by_name(dynamic_cfg, :l_spring, type)
    p0 = get_ring_prop_by_name(dynamic_cfg, :p0, type)
    get_area0(num_particles, l_spring, p0)
end
get_area0(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_area0)

"""
Returns the p0 at which `area0` is equal to the equilibrium area
considering only the springs.
"""
function get_equilibrium_p0(num_particles::Int) 
    theta = 2 * π  / num_particles
    return 2 * (num_particles * (1 - cos(theta))/sin(theta))^.5
end
get_equilibrium_p0(dynamic_cfg::RingsCfg{U, T, I}, type=nothing) where {U, T, I} = get_equilibrium_p0(ring_num_particles(dynamic_cfg, type))
get_equilibrium_p0(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_equilibrium_p0)


"""
Equilibrium area of the polygon formed by the centers of the ring particles.

NOTE: The area contribution from the particles is not considered here.
"""
function get_equilibrium_area(dynamic_cfg::RingsCfg{U, T, I}, type=nothing) where {U, T, I}
    l_spring = get_ring_prop_by_name(dynamic_cfg, :l_spring, type)
    num_particles = ring_num_particles(dynamic_cfg, type) 
    a0 = get_area0(dynamic_cfg, type)

    p0_lim = get_equilibrium_p0(dynamic_cfg, type)
    a0_lim = (num_particles * l_spring / p0_lim)^2 

    if a0 < a0_lim
        return a0
    end
    theta = 2 * π / num_particles

    k_a = get_ring_prop_by_name(dynamic_cfg, :k_area, type)
    k_m = get_ring_prop_by_name(dynamic_cfg, :k_spring, type)
    angle = π * (1 - 2 / num_particles)

    get_area(l) = num_particles * l^2 / (4 * tan(π / num_particles))
    function get_fa(l)
        a = get_area(l)
        return k_a * (a0 - a) * l * sin(angle / 2)
    end

    get_fl(l) = k_m * (l - l_spring)

    function func!(F, l_sqrt)
        l = l_sqrt[1]^2
        F[1] = get_fa(l) -  2 * get_fl(l)
    end
    sol = nlsolve(func!, [sqrt(l_spring * 1.1)])
    l_sol = sol.zero[1]^2
    return get_area(l_sol)

    # function get_r(f)
    #     sqrt(f * 2 * a0 / (num_particles * sin(theta)))
    # end

    # function get_fm(f)
    #     k_m * (get_r(f) * sqrt(2 * (1 - cos(theta))) - l_spring)
    # end

    # function get_fm_total(f)
    #     2 * get_fm(f) * sin(theta / 2)
    # end

    # function get_fa(f)
    #     r = get_r(f)
    #     k_a * (a0 - 10/2 * r^2 * sin(theta)) * r * sin(theta)
    # end

    # function func!(F, f)
    #     F[1] = get_fa(f[1]^2) - get_fm_total(f[1]^2)
    # end

    # sol = nlsolve(func!, [sqrt(0.5)])
    # f_sol = sol.zero[1]^2
    # a0_sol = f_sol * a0 
    # return a0_sol
end
get_equilibrium_area(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_equilibrium_area)

"""
Contribution of the particles to the area of the rings.
The area of a ring is the area of the polygon formed by the centers of its particles (A_p),
plus the area of the particles that lies outside this polygon (A_c); this
function returns A_c.
"""
function get_particles_area_contribution(dynamic_cfg::RingsCfg{U, T, I}, type=nothing) where {U, T, I}
    n = ring_num_particles(dynamic_cfg, type)
    l_spring = get_ring_prop_by_name(dynamic_cfg, :l_spring, type)
    diameter = get_particle_radius(dynamic_cfg, type) * 2

    root_term = (diameter^2 - l_spring^2)^.5
    t = π/2 - atan(l_spring/root_term)
    area_intersect = 1/4 * (diameter^2 * t - l_spring * root_term)

    return n * π / 4 * diameter^2 * (1 - (n-2)/(2*n)) - n * area_intersect
end
get_particles_area_contribution(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_particles_area_contribution)

"""
Equilibrium ring radius, calculated as follow: assuming the ring has its equilibrium area and 
its shape is a regular polygon, the radius is the segment length from the center to a vertex plus the particle radius.
"""
function get_ring_radius(dynamic_cfg::RingsCfg{U, T, I}, type=nothing) where {U, T, I}
    area_eq = get_equilibrium_area(dynamic_cfg, type)
    p_radius = get_particle_radius(dynamic_cfg, type)
    num_particles = ring_num_particles(dynamic_cfg, type)
    
    radius_to_particle = (2 * area_eq / (num_particles * sin(2 * pi / num_particles)))^.5
    return radius_to_particle + p_radius
end
get_ring_radius(dynamic_cfg::RingsCfg{U, T, I}) where {U<:AbstractVector, T, I} = get_rings_property(dynamic_cfg, get_ring_radius)

# =
# Integration Configs
# = 

struct InvasionsCfg
    steps_to_update::Int
end

struct IntCfgExtra{RC<:Union{ChunksCfg, Nothing}, InvT<:Union{InvasionsCfg, Nothing}}
    r_chunks_cfg::RC
    invasions_cfg::InvT
end
StructTypes.StructType(::Type{D}) where D <: IntCfgExtra = StructTypes.Struct()

function RingsIntCfg(; dt, p_chunks_cfg=nothing, r_chunks_cfg=nothing, 
    invasions_cfg=nothing, device=nothing) 
    if isnothing(device)
        device = MaviCfg.IntCfg(dt=dt).device
    end

    MaviCfg.IntCfg(dt, p_chunks_cfg, device, IntCfgExtra(r_chunks_cfg, invasions_cfg))
end
function RingsIntCfg(int_cfg::MaviCfg.IntCfg; r_chunks_cfg=nothing, 
    invasions_cfg=nothing) 
    RingsIntCfg(
        dt=int_cfg.dt,
        p_chunks_cfg=int_cfg.chunks_cfg,
        device=int_cfg.device,
        r_chunks_cfg=r_chunks_cfg,
        invasions_cfg=invasions_cfg,
    )
end

end