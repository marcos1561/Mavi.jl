module RingStateTests

using Test, Random
using Mavi
using Mavi.StaticArrays
using Mavi.Rings: RingsState, ActiveState
using Mavi.Rings.Configs: RingsCfg
using Mavi.Rings.States: add_ring!, remove_ring!, get_rings_ids, get_ring_particles_indices
using Mavi.States: get_particles_ids
using Mavi.Systems: is_valid_pair, get_entity_vel

function check_valid_pairs()
    pos = fill(SVector(0.0, 0.0), 4, 3)
    states = (
        RingsState(rings_pos=copy(pos), pol=zeros(3)),
        RingsState(rings_pos=copy(pos), pol=zeros(3),
            active_state=ActiveState([true, false, true])),
        RingsState(rings_pos=copy(pos), pol=zeros(3),
            types=[1, 2, 1], num_particles=[3, 4]),
        RingsState(rings_pos=copy(pos), pol=zeros(3),
            types=[1, 2, 1], num_particles=[3, 4],
            active_state=ActiveState([true, false, true])),
        RingsState(rings_pos=copy(pos), pol=zeros(3),
            active_state=ActiveState(fill(false, 3))),
        RingsState(rings_pos=Matrix{SVector{2, Float64}}(undef, 4, 0), pol=Float64[]),
    )

    for state in states
        valid_ids = Set(get_particles_ids(state))
        system = (state=state, dynamic_cfg=nothing)
        for i in (-1, 0, 1, 3, 4, 5, 8, 9, 11, 12, 13),
            j in (-1, 0, 1, 3, 4, 5, 8, 9, 11, 12, 13)
            expected = i in valid_ids && j in valid_ids
            @test is_valid_pair(state, nothing, i, j) == expected
            @test is_valid_pair(state, i, j) == expected
            @test is_valid_pair(system, i, j) == expected
        end
    end

    # The mask changes immediately, even when ID lists are rebuilt in a batch.
    state = RingsState(rings_pos=copy(pos), pol=zeros(3),
        active_state=ActiveState([true, false, true]))
    remove_ring!(state, 1)
    @test !is_valid_pair(state, nothing, 1, 9)
    @test !is_valid_pair(state, nothing, 9, 1)
    @test is_valid_pair(state, nothing, 9, 10)
    @test add_ring!(state, pos[:, 1], 0.0) == 1
    @test is_valid_pair(state, nothing, 1, 9)
    @test add_ring!(state, pos[:, 2], 0.0) == 2
    @test is_valid_pair(state, nothing, 5, 9)
end

function entity_velocity_allocations(system, id)
    get_entity_vel(system, id)
    return @allocated get_entity_vel(system, id)
end

function check_entity_velocities()
    for T in (Float32, Float64), typed in (false, true), variable in (false, true)
        @testset "$T, typed=$typed, variable=$variable" begin
            num_particles = typed ? [3, 5] : 3
            num_max = typed ? 5 : 3
            state = RingsState(
                rings_pos=fill(zero(SVector{2, T}), num_max, 4),
                pol=T[0, pi/2, pi, -pi/2],
                types=typed ? [2, 1, 2, 1] : nothing,
                num_particles=num_particles,
                active_state=variable ? ActiveState([true, false, true, true]) : nothing,
            )
            potential = Mavi.Configs.HarmTruncCfg(
                k_rep=T(1), k_atr=T(0), dist_eq=T(1), dist_max=T(1.2),
            )
            cfg = RingsCfg(
                p0=T(3.5), relax_time=T(1), rot_diff=T(0),
                k_area=T(1), k_spring=T(1), l_spring=T(1),
                vo=typed ? T[1, 2] : T(2),
                mobility=typed ? T[0.5, 1.5] : T(1.5),
                num_particles=num_particles, interaction_finder=potential,
            )
            # Padding and inactive rings must not contribute to the velocity.
            forces = fill(SVector{2, T}(NaN, NaN), length(state.pos))
            for id in get_rings_ids(state)
                particle_ids = get_ring_particles_indices(state, id)
                for pid in particle_ids
                    forces[pid] = zero(SVector{2, T})
                end
                forces[first(particle_ids)] = length(particle_ids) * SVector{2, T}(2, -1)
            end
            system = (
                state=state, dynamic_cfg=cfg, forces=(forces,),
                int_cfg=(dt=T(0.01),), rng=Xoshiro(123),
            )
            expected = typed ?
                [SVector(5, -1.5), SVector(1, 0.5), SVector(1, -1.5), SVector(1, -1.5)] :
                [SVector(5, -1.5), SVector(3, 0.5), SVector(1, -1.5), SVector(3, -3.5)]

            for id in get_rings_ids(state)
                velocity = @inferred get_entity_vel(system, id)
                @test velocity isa SVector{2, T}
                @test velocity ≈ SVector{2, T}(expected[id])
                @test entity_velocity_allocations(system, id) == 0
            end
            @test_throws BoundsError get_entity_vel(system, 0)
            @test_throws BoundsError get_entity_vel(system, 5)

            # Predict the center displacement of the actual integration kernel.
            old_pos = copy(state.pos)
            Mavi.Rings.Integration.update!(system)
            for id in get_rings_ids(state)
                particle_ids = get_ring_particles_indices(state, id)
                displacement = sum(state.pos[pid] - old_pos[pid] for pid in particle_ids)
                measured_velocity = displacement / (length(particle_ids) * system.int_cfg.dt)
                @test measured_velocity ≈ SVector{2, T}(expected[id])
            end
        end
    end
end

end
