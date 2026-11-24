module SourceTests

using Test
using Mavi
using Mavi.StaticArrays
using Mavi.ChunksMod: Chunks
using Mavi.Configs: RectangleCfg, SpaceCfg, PotentialWalls, HarmTruncCfg
using Mavi.Rings.Sources: ChunksChecker

function check_boundary_chunks()
    cases = [
        (offset=(-0.5, -0.5), size=(1.0, 1.0), rows=4:4, cols=1:1),
        (offset=(7.5, 7.5), size=(1.0, 1.0), rows=1:1, cols=4:4),
        (offset=(-1.0, -1.0), size=(10.0, 10.0), rows=1:4, cols=1:4),
        (offset=(1.5, 1.5), size=(1.0, 1.0), rows=3:4, cols=1:2),
        (offset=(2.0, 2.0), size=(0.0, 0.0), rows=3:4, cols=1:2),
        (offset=(-2.0, -2.0), size=(1.0, 1.0), rows=1:0, cols=1:0),
        (offset=(9.0, 9.0), size=(1.0, 1.0), rows=1:0, cols=1:0),
        (offset=(0.0, 0.0), size=(8.0, 8.0), rows=1:4, cols=1:4),
    ]

    for origin in (SVector(0.0, 0.0), SVector(3.0, -2.0))
        space = SpaceCfg(
            geometry_cfg=RectangleCfg(length=8.0, height=8.0, bottom_left=origin),
            wall_type=PotentialWalls(potential=HarmTruncCfg(
                k_rep=30.0, k_atr=0.0, dist_eq=0.5, dist_max=0.5,
            )),
        )
        pos = SVector{2, Float64}[]
        chunks = Chunks(4, 4, space, pos, 0.5)
        boxes = [RectangleCfg(
            length=case.size[1], height=case.size[2],
            bottom_left=origin + SVector(case.offset),
        ) for case in cases]
        checker = ChunksChecker(pos, chunks, boxes)

        for (case, ids) in zip(cases, checker.ids)
            expected = Set(CartesianIndex(row, col) for col in case.cols for row in case.rows)
            @test Set(ids) == expected
            @test all(id -> checkbounds(Bool, chunks.num_particles_in_chunk, id), ids)
        end
    end
end

end
