"""
Example: Rings with a sink and a source

Creates an empty system with a source and a sink at apposite sides.
"""
module Example

using Mavi.StaticArrays

using Mavi.Rings
using Mavi.Rings.Configs
using Mavi.Rings.InitStates
using Mavi.Rings.Sources

using Mavi.ChunksMod
using Mavi.Visualization
using Mavi.Visualization.SystemGraphs

struct SinkSourceGraphCfg <: GraphCompCfg end
struct SinkSourceGraph <: GraphCompDebug end
function SystemGraphs.get_graph(ax, pos_obs, system, cfg::SinkSourceGraphCfg)
    source = system.info.sources[1]
    sink = system.info.sources[2]

    for idx in eachindex(source.bbox)
        bbox = source.bbox[idx]
        drawn_borders(ax, bbox, adjust_lims=false, color=:green)
    end

    drawn_borders(ax, sink.cfg.geometry_cfg, adjust_lims=false, color=:red)

    return SinkSourceGraph()
end

function create_system(;num_particles=10)
    interaction_cfg = Configs.HarmTruncCfg(
        k_rep = 30.,
        k_atr = 3,
        dist_eq = 1.,
        dist_max = 1 * 1.5,
    )

    dynamic_cfg = Configs.RingsCfg(
        p0=3.545,
        relax_time=100.,
        vo=1.,
        mobility=1.,
        rot_diff=0.,
        k_area=3.,
        k_spring=30.,
        num_particles=num_particles,
        l_spring=interaction_cfg.dist_eq*0.8,
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
        pol=random_pol(num_rings),
        active_state=ActiveState(fill(false, num_rings)),
    )

    bbox = Configs.get_bounding_box(space_cfg.geometry_cfg)
    max_size = interaction_cfg.dist_max * 1.2
    num_x_chunks = floor(Int, bbox.length / max_size)
    num_y_chunks = floor(Int, bbox.height / max_size)
    
    system = RingsSystem(
        state = state,
        space_cfg = space_cfg,
        dynamic_cfg = dynamic_cfg,
        int_cfg = Configs.RingsIntCfg(
            dt=0.01,
            p_chunks_cfg=Configs.ChunksCfg(
                num_x_chunks, num_y_chunks,
            ),
            device=Configs.Threaded(),
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
    )

    return system
end

function main(test=false)
    system = create_system()
    
    anim_cfg = AnimationCfg(
        num_steps_per_frame=30,
        graph_cfg=MainGraphCfg([
            CircleGraphCfg(),
            SinkSourceGraphCfg(),
        ])
    )

    if !test
        animate(system, anim_cfg)
    else
        Rings.Mavi.run(system, tf=1)
    end
end

end

if !((@isdefined TEST_EX) && TEST_EX)
    Example.main()
end