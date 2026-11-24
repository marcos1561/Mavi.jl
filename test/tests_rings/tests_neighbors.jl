module NeighborTests

using Test, Random
using Mavi
using Mavi.StaticArrays
using Mavi.Configs: Sequencial, Threaded, RectangleCfg, SpaceCfg, PotentialWalls, ChunksCfg
using Mavi.Rings: RingsSystem, RingsState, ActiveState
using Mavi.Rings.Configs: RingsCfg, RingsIntCfg, HarmTruncCfg
using Mavi.Rings.States: get_rings_ids, get_ring_particles_indices, remove_ring!
using Mavi.Rings.NeighborsMod
using Mavi.MaviSerder: JSON3

function record_on_workers(f, device)
    if device isa Threaded
        Threads.@threads :static for _ in 1:(4 * Threads.nthreads())
            f()
        end
    else
        for _ in 1:4
            f()
        end
    end
end

function check_repeatable_neighbors()
    for device in (Sequencial(), Threaded()), only_count in (false, true)
        @testset "$(typeof(device)), only_count=$only_count" begin
            cfg = NeighborsCfg(only_count=only_count, tol=1.25)
            neigh = RepeatableNeighbors(
                num_entities=18, num_max_neighbors=only_count ? nothing : 1,
                device=device, cfg=cfg,
            )
            # Every worker observes the same pairs, in both orientations.
            record_on_workers(device) do
                for j in 2:17
                    neigh_update!(neigh, 1, j, 1.1, 1.0)
                    neigh_update!(neigh, j, 1, 1.1, 1.0)
                end
                neigh_update!(neigh, 1, 1, 0.0, 1.0)
                neigh_update!(neigh, 1, 18, 1.25, 1.0)
            end
            neigh_sum_buffers(neigh)
            expected = [16; ones(Int, 16); 0]
            @test get_neigh_count(neigh) == expected
            if !only_count
                @test get_neigh_list(neigh, 1) == collect(2:17)
                @test all(i -> get_neigh_list(neigh, i) == [1], 2:17)
                @test isempty(get_neigh_list(neigh, 18))
                @test size(get_neigh(neigh).list, 1) >= 16
            else
                @test isnothing(get_neigh(neigh).list)
            end
            @test_throws BoundsError neigh_update!(neigh, 1, 19, 0.5, 1.0)
            neigh_sum_buffers(neigh)
            @test get_neigh_count(neigh) == expected

            neigh_clean!(neigh)
            @test all(iszero, get_neigh_count(neigh))
            record_on_workers(device) do
                neigh_update!(neigh, 2, 3, 0.5, 1.0)
            end
            neigh_sum_buffers(neigh)
            @test get_neigh_count(neigh) == [0; 1; 1; zeros(Int, 15)]
            if !only_count
                @test isempty(get_neigh_list(neigh, 1))
                @test get_neigh_list(neigh, 2) == [3]
                @test get_neigh_list(neigh, 3) == [2]
            end
            neigh_clean!(neigh)
            neigh_sum_buffers(neigh)
            @test all(iszero, get_neigh_count(neigh))
        end
    end

    @testset "Empty and ordinary neighbors" begin
        for only_count in (false, true)
            neigh = RepeatableNeighbors(num_entities=0,
                num_max_neighbors=only_count ? nothing : 0,
                device=Sequencial(), cfg=NeighborsCfg(only_count=only_count))
            neigh_sum_buffers(neigh)
            @test isempty(get_neigh_count(neigh))
        end
        ordinary = Neighbors(num_entities=2, num_max_neighbors=2,
            device=Sequencial(), cfg=NeighborsCfg())
        for _ in 1:2
            neigh_update!(ordinary, 1, 2, 0.5, 1.0)
        end
        neigh_sum_buffers(ordinary)
        @test get_neigh_count(ordinary) == [2, 2]
        @test get_neigh_list(ordinary, 1) == [2, 2]
    end
end

function make_system(device, only_count; use_chunks=true)
    centers = [SVector(3.0, 3.0), SVector(3.4, 3.4), SVector(4.0, 3.2), SVector(12.0, 12.0)]
    pos = [center + 0.3 * SVector(cos((i-1)*pi/2), sin((i-1)*pi/2))
        for i in 1:4, center in centers]
    state = RingsState(rings_pos=pos, pol=zeros(4),
        active_state=ActiveState([true, false, true, true]))
    potential = HarmTruncCfg(k_rep=30.0, k_atr=3.0, dist_eq=1.0, dist_max=1.5)
    RingsSystem(
        state=state,
        space_cfg=SpaceCfg(geometry_cfg=RectangleCfg(length=20.0, height=20.0),
            wall_type=PotentialWalls(potential=potential)),
        dynamic_cfg=RingsCfg(p0=3.5, relax_time=1.0, vo=1.0, mobility=1.0,
            rot_diff=0.0, k_area=1.0, k_spring=2.0, l_spring=0.8,
            num_particles=4, interaction_finder=potential),
        int_cfg=RingsIntCfg(dt=0.01, device=device,
            p_chunks_cfg=use_chunks ? ChunksCfg(10, 10) : nothing),
        r_neighbors_cfg=NeighborsCfg(only_count=only_count),
        rng=MersenneTwister(1234),
    )
end

function check_graph(system)
    neigh = system.info.r_neigh
    @test neigh isa RepeatableNeighbors
    expected = [Int[] for _ in axes(system.state.rings_pos, 2)]
    ids = collect(get_rings_ids(system.state))
    cutoff = 1.0 * get_neigh(neigh).cfg.tol
    for (k, ri) in enumerate(ids), rj in ids[k+1:end]
        pi = get_ring_particles_indices(system.state, ri)
        pj = get_ring_particles_indices(system.state, rj)
        if any(sum(abs2, system.state.pos[i] - system.state.pos[j]) < cutoff^2 for i in pi, j in pj)
            push!(expected[ri], rj)
            push!(expected[rj], ri)
        end
    end
    @test get_neigh_count(neigh) == length.(expected)
    if !get_neigh(neigh).cfg.only_count
        for i in eachindex(expected)
            @test get_neigh_list(neigh, i) == sort(expected[i])
        end
    end
end

function check_ring_neighbors()
    for (device, use_chunks) in ((Sequencial(), false), (Sequencial(), true), (Threaded(), true)),
        only_count in (false, true)
        system = make_system(device, only_count; use_chunks=use_chunks)
        check_graph(system)
        @test get_neigh_count(system.info.r_neigh) == [1, 0, 1, 0]
        for _ in 1:3
            Mavi.Rings.Integration.step!(system)
            check_graph(system)
        end
        copied = deepcopy(system)
        check_graph(copied)
        @test get_neigh(copied.info.r_neigh).cfg.only_count == only_count
        remove_ring!(system.state, 3)
        Mavi.States.update_ids!(system)
        Mavi.Rings.Integration.system_initialization(system)
        check_graph(system)
        @test all(iszero, get_neigh_count(system.info.r_neigh))
    end

    @testset "Serialization" begin
        for only_count in (false, true)
            system = make_system(Sequencial(), only_count)
            mktempdir() do root
                Mavi.save_system(system, root)
                restored = Mavi.load_system(root)
                check_graph(restored)
                @test get_neigh(restored.info.r_neigh).cfg.only_count == only_count
                Mavi.Rings.Integration.step!(restored)
                check_graph(restored)

                # Older snapshots name RingsInfo with a plain Neighbors field.
                info = system.info
                legacy_info = Mavi.Rings.RingsInfo(info.continuos_pos, info.areas,
                    info.cms, info.p_neigh, get_neigh(info.r_neigh), info.r_chunks,
                    info.invasions, info.sources, info.user_data)
                path = joinpath(root, "configs.json")
                data = JSON3.read(read(path, String), Dict{String, Any})
                data["info"]["type"] = string(typeof(legacy_info))
                open(path, "w") do io
                    JSON3.write(io, data)
                end
                check_graph(Mavi.load_system(root))
            end
        end
    end
end

end
