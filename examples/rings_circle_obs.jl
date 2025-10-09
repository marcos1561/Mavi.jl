"""
Example: Rings with circular obstacles

Rings are created in a rectangular grid, then the ones intersecting with obstacles
are removed. Two circular obstacles, with slippery walls, are created. 
"""
module Example

using Mavi.Rings
using Mavi.Visualization

"Filter out every ring inside `remove_region`."
function filter_rings_pos(rings_pos, remove_region)
    filtered_ids = Int[]
    remove_region_bbox = Configs.get_bounding_box(remove_region)
    for ring_id in axes(rings_pos, 2)
        ring_bbox = Configs.RectangleCfg(rings_pos[:, ring_id])

        if !Configs.check_intersection(ring_bbox, remove_region_bbox)
            push!(filtered_ids, ring_id)
        end
    end
    
    return rings_pos[:, filtered_ids]
end

function create_system(;num_particles, num_cols, num_rows)
    interaction_cfg = Configs.HarmTruncCfg(
        k_rep=40,
        k_atr=4,
        dist_eq=1,
        dist_max=1 + 0.2,
    )

    dynamic_cfg = Configs.RingsCfg(
        p0=3.5,
        relax_time=1.0,
        vo=1.0,
        mobility=1.0,
        rot_diff=0.05,
        k_area=1.0,
        k_spring=40.0,
        l_spring=1.0,
        num_particles=num_particles,
        interaction_finder=interaction_cfg,
    )

    num_rings = num_cols * num_rows
    rings_pos, geometry_cfg = InitStates.rectangular_grid(
        num_cols=num_cols,
        num_rows=num_rows,
        num_particles=num_particles,
        p_radius=Configs.particle_radius(dynamic_cfg),
        pad_x=0.1, pad_y=0.1,
        to_svector=true,
    )

    ring_r = Configs.get_ring_radius(dynamic_cfg)
    bl = geometry_cfg.bottom_left
    length, height = geometry_cfg.size

    circle1_cfg = Configs.CircleCfg(
        radius=2*ring_r,
        center=[bl.x + length/2, bl.y + 1/4 * height],
    )

    circle2_cfg = Configs.CircleCfg(
        radius=2*ring_r,
        center=[bl.x + length/2, bl.y + (1 - 1/4) * height],
    )

    rings_pos = filter_rings_pos(rings_pos, circle1_cfg)
    rings_pos = filter_rings_pos(rings_pos, circle2_cfg)

    space_cfg = Configs.SpaceCfg([
        (Configs.PeriodicWalls(), geometry_cfg),
        (Configs.SlipperyWalls(), circle1_cfg),
        (Configs.SlipperyWalls(), circle2_cfg),
    ])

    state = RingsState(
        rings_pos=rings_pos,
        pol=InitStates.random_pol(num_rings),
    )

    bbox = space_cfg.geometry_cfg.list[1]
    # bbox = space_cfg.geometry_cfg
    
    max_size = interaction_cfg.dist_max * 1.1
    num_x_chunks = floor(Int, bbox.length / max_size)
    num_y_chunks = floor(Int, bbox.height / max_size)

    system = RingsSystem(
        state=state,
        space_cfg=space_cfg,
        dynamic_cfg=dynamic_cfg,
        int_cfg=Configs.RingsIntCfg(
            dt=0.01,
            p_chunks_cfg=Configs.ChunksCfg(num_x_chunks, num_y_chunks),
            # device=Configs.Threaded(),
        ),
    )

    return system
end

function main(test=false)
    system = create_system(
        num_particles=10,
        num_cols=13,
        num_rows=13,
    )

    anim_cfg = AnimationCfg(
        num_steps_per_frame=15,
        graph_cfg=CircleGraphCfg(),
        # begin_paused=true,
        # graph_cfg=ScatterGraphCfg(),
    )
    
    if !test
        animate(system, anim_cfg)
    else
        Rings.Mavi.run_system(system, tf=1)
    end
end

end

if !((@isdefined TEST_EX) && TEST_EX)
    Example.main()
end