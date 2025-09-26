module RingsUtils

export state_square_distance
export create_system_normal, create_system_source, create_system_types
export create_systems_list

using StaticArrays, Random

using Mavi.Rings
using Mavi.Rings.Sources
using Mavi.Rings.States
using Mavi.Rings.Configs
using Mavi.Rings.InitStates

function state_square_distance(s1, s2)
    dist = 0
    for p_id in get_particles_ids(s1)
        dist += sum((s1.state.pos[p_id] - s2.state.pos[p_id]).^2)
    end
    for r_id in get_rings_ids(s1)
        dist += sum((s1.state.pos[r_id] - s2.state.pos[r_id]).^2)
    end

    return dist
end

function create_system_normal(rng=nothing, device=Configs.Sequencial(), use_chunks=true)
    num_particles = 10
    num_cols = 13 
    num_rows = 13

    interaction_cfg = HarmTruncCfg(
        k_rep=20,
        k_atr=4,
        dist_eq=1,
        dist_max=1 + 0.2,
    )

    dynamic_cfg = RingsCfg(
        p0=3.5,
        relax_time=1,
        vo=1.0,
        mobility=1,
        rot_diff=0.05,
        k_area=1,
        k_spring=20,
        l_spring=1,
        num_particles=num_particles,
        interaction_finder=interaction_cfg,
    )

    num_rings = num_cols * num_rows
    rings_pos, geometry_cfg = rectangular_grid(
        num_cols=num_cols,
        num_rows=num_rows,
        num_particles=num_particles,
        p_radius=Configs.particle_radius(dynamic_cfg),
        pad_x=0.1, pad_y=0.1,
    )

    space_cfg = Configs.SpaceCfg(
        geometry_cfg=geometry_cfg,
        wall_type=Configs.PeriodicWalls(),
    )
    
    state = RingsState(
        rings_pos=rings_pos,
        pol=random_pol(num_rings, rng=rng),
    )

    bbox = space_cfg.geometry_cfg
    max_size = interaction_cfg.dist_max * 1.1
    num_x_chunks = floor(Int, bbox.length / max_size)
    num_y_chunks = floor(Int, bbox.height / max_size)

    p_chunks_cfg = nothing
    if use_chunks
        p_chunks_cfg=Configs.ChunksCfg(num_x_chunks, num_y_chunks)
    end

    system = RingsSystem(
        state=state,
        space_cfg=space_cfg,
        dynamic_cfg=dynamic_cfg,
        int_cfg=RingsIntCfg(
            dt=0.01,
            p_chunks_cfg=p_chunks_cfg,
            device=device,
        ),
        rng=rng
    )

    return system
end

function create_system_source(rng=nothing, device=Configs.Sequencial(), use_chunks=true)
    num_particles=10

    interaction_cfg = Configs.HarmTruncCfg(
        k_rep = 30,
        k_atr = 3,
        dist_eq = 1,
        dist_max = 1 * 1.5,
    )

    dynamic_cfg = Configs.RingsCfg(
        p0=3.545,
        relax_time=100,
        vo=1,
        mobility=1,
        rot_diff=0,
        k_area=3,
        k_spring=30,
        l_spring=interaction_cfg.dist_eq*0.8,
        num_particles=num_particles,
        interaction_finder=interaction_cfg,
    )

    ring_d = get_ring_radius(dynamic_cfg) * 2

    geometry_cfg = Configs.RectangleCfg(
        length=18 * ring_d,
        height=10 * ring_d,
    )

    rings_pos = Array{Float64, 3}(undef, 2, num_particles, 100)

    ring_area = get_equilibrium_area(dynamic_cfg)
    r = (ring_area / π)^.5
    spawn_pos = create_circle((0, 0), r, num_particles)
    spawn_pos = convert(Matrix{Float64}, spawn_pos)
    spawn_pos = copy(reinterpret(SVector{2, Float64}, vec(spawn_pos)))

    space_cfg = Configs.SpaceCfg(
        geometry_cfg = geometry_cfg,
        wall_type = Configs.PeriodicWalls(),
    )
        
    num_rings = size(rings_pos, 3)
    state = RingsState(
        rings_pos=rings_pos,
        pol=random_pol(num_rings, rng=rng),
        active_state=ActiveState(fill(false, num_rings)),
    )

    bbox = Configs.get_bounding_box(space_cfg.geometry_cfg)
    max_size = interaction_cfg.dist_max * 1.2
    num_x_chunks = floor(Int, bbox.length / max_size)
    num_y_chunks = floor(Int, bbox.height / max_size)
    
    p_chunks_cfg = nothing
    if use_chunks
        p_chunks_cfg=Configs.ChunksCfg(
            num_x_chunks, num_y_chunks,
        )
    end

    system = RingsSystem(
        state = state,
        space_cfg = space_cfg,
        dynamic_cfg = dynamic_cfg,
        int_cfg = RingsIntCfg(
            dt=0.01,
            p_chunks_cfg=p_chunks_cfg,
            device=device,
        ),
        source_cfg=[
            SourceCfg(
                bottom_left=geometry_cfg.bottom_left + SVector(2*ring_d, geometry_cfg.height/2 - ring_d),
                spawn_pos=spawn_pos,
                size=(2, 2),
                spawn_pol=:random,
                pad=get_particle_radius(dynamic_cfg),
            ),
            SinkCfg(Configs.RectangleCfg(
                length=ring_d*2,
                height=ring_d*2,
                bottom_left=geometry_cfg.bottom_left + SVector(geometry_cfg.length - 4*ring_d, geometry_cfg.height/2 - ring_d),
            )),
        ],
        rng=rng,
    )

    return system
end

function create_system_types(rng=nothing, device=Configs.Sequencial(), use_chunks=true)
    if isnothing(rng)
        rng = Random.GLOBAL_RNG
    end

    num_cols = 5
    num_rows = 5
    interaction_cfg_1 = HarmTruncCfg(
        k_rep=13,
        k_atr=0.1,
        dist_eq=1,
        dist_max=1 * 1.1,
    )

    interaction_cfg_2 = HarmTruncCfg(
        k_rep=13,
        k_atr=0.1,
        dist_eq=0.8,
        dist_max=0.8 * 1.1,
    )

    pair_dist_eq = interaction_cfg_1.dist_eq/2 + interaction_cfg_2.dist_eq/2
    interaction_cfg_pair = HarmTruncCfg(
        k_rep = 13,
        k_atr = 10,
        dist_eq = pair_dist_eq,
        dist_max = pair_dist_eq * 1.2,
    )
    
    interaction_finder = InteractionMatrix([
        [interaction_cfg_1, interaction_cfg_pair];; 
        [interaction_cfg_pair, interaction_cfg_2]
    ])

    interactions = list_interactions(interaction_finder)
    self_interactions = list_self_interactions(interaction_finder)

    num_particles = [10, 5]

    dynamic_cfg = RingsCfg(
        p0=3.5,
        relax_time=1,
        vo=1,
        mobility=1,
        rot_diff=0.5,
        k_area=1,
        k_spring=20,
        l_spring=[inter.dist_eq*0.8 for inter in self_interactions],
        num_particles=num_particles,
        interaction_finder=interaction_finder,
    )

    interactions = list_interactions(dynamic_cfg.interaction_finder)

    num_rings = num_cols * num_rows
    types = rand(rng, 1:2, num_rings)
    p_radius = Configs.particle_radius.([interaction_cfg_1, interaction_cfg_2])
    rings_pos, geometry_cfg = InitStates.rectangular_grid(
        num_cols=num_cols,
        num_rows=num_rows,
        num_particles=dynamic_cfg.num_particles,
        p_radius=p_radius,
        types=types,
        pad_x=0.1, pad_y=0.1,
    )

    space_cfg = Configs.SpaceCfg(
        geometry_cfg=geometry_cfg,
        wall_type=Configs.PeriodicWalls(),
    )
    
    state = RingsState(
        rings_pos=rings_pos,
        pol=InitStates.random_pol(num_rings, rng=rng),
        types=types,
        num_particles=num_particles
    )

    bbox = Configs.get_bounding_box(space_cfg.geometry_cfg)
    max_size = maximum([i.dist_max for i in interactions]) * 1.1
    num_x_chunks = floor(Int, bbox.length / max_size)
    num_y_chunks = floor(Int, bbox.height / max_size)

    p_chunks_cfg = nothing
    if use_chunks
        p_chunks_cfg=Configs.ChunksCfg(
            num_x_chunks, num_y_chunks,
        )
    end

    system = RingsSystem(
        state=state,
        space_cfg=space_cfg,
        dynamic_cfg=dynamic_cfg,
        int_cfg=RingsIntCfg(
            dt=0.01,
            p_chunks_cfg=p_chunks_cfg,
            device=device,
        ),
        rng=rng,
    )

    return system
end

create_systems_list = (
    normal=create_system_normal,
    source=create_system_source,
    types=create_system_types,
)

end # RingsUtils

