module RingStateTests

using Test
using Mavi
using Mavi.StaticArrays
using Mavi.Rings: RingsState, ActiveState
using Mavi.Rings.States: add_ring!, remove_ring!
using Mavi.States: get_particles_ids
using Mavi.Systems: is_valid_pair

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

end
