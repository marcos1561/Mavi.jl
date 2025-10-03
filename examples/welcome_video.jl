"""
Example: Welcome Video
"""
module Example

using Mavi.StaticArrays

using Mavi.Rings
using Mavi.Rings.Configs
using Mavi.Rings.Integration
using Mavi.Rings.InitStates
using Mavi.Rings.Sources

using Mavi.ChunksMod
using Mavi.Visualization
using Mavi.Visualization.SystemGraphs

struct SinkSourceGraphCfg <: GraphCompCfg end
struct SinkSourceGraph <: GraphCompDebug end
"Show sink source borders."
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
        rot_diff=0.1,
        k_area=3.,
        k_spring=30.,
        l_spring=interaction_cfg.dist_eq*0.8,
        num_particles=num_particles,
        interaction_finder=interaction_cfg,
    )

    ring_d = get_ring_radius(dynamic_cfg) * 2

    # =
    # Space Cfg
    # =

    # Mavi letters 
    lines_list = [] 
    circles_list = []
    
    # Pad from borders
    mavi_pad = SVector(4 * ring_d, 2 * ring_d)
    
    # Letter configuration
    w = 2*ring_d # Letter A width (The letter M has double this width)  
    h = 4*ring_d # Letter hight
    pad = w/1.5 # Horizontal pad between letters
    
    circle_r = w / 8

    # Letter base point
    bp = mavi_pad

    # M
    push!(lines_list, [bp, bp + SVector(w/2, h)])
    push!(lines_list, [bp + SVector(w/2, h), bp + SVector(w, 0)])
    push!(lines_list, [bp + SVector(w, 0), bp + SVector(w/2 + w, h)])
    push!(lines_list, [bp + SVector(w/2 + w, h), bp + SVector(2*w, 0)])
    
    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp,
    ))
    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp + SVector(2*w, 0),
    ))

    bp = bp + SVector(2*w + pad, 0)
    
    # A
    push!(lines_list, [bp, bp + SVector(w/2, h)])
    push!(lines_list, [bp + SVector(w/2, h), bp + SVector(w, 0)])
    push!(lines_list, [bp + SVector(w/4, 2*h * 1/4), bp + SVector(3/4*w, 2*h * 1/4)])

    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp,
    ))
    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp + SVector(w, 0),
    ))
    
    bp = bp + SVector(w + pad, 0)
    
    # V
    push!(lines_list, [bp + SVector(0, h), bp + SVector(w/2, 0)])
    push!(lines_list, [bp + SVector(w, h), bp + SVector(w/2, 0)])
    
    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp + SVector(0, h),
    ))
    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp + SVector(w, h),
    ))
    
    
    bp = bp + SVector(w + pad, 0)
    
    # I
    push!(lines_list, [bp, bp + SVector(0, h)])
    
    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp,
    ))
    push!(circles_list, Configs.CircleCfg(
        radius=circle_r,
        center=bp + SVector(0, h),
    ))
    

    lines_cfg = Configs.LinesCfg(lines_list)

    # Rectangle border
    mavi_bbox = Configs.get_bounding_box(lines_cfg)

    rect_cfg = Configs.RectangleCfg(
        length=2*mavi_pad[1] + mavi_bbox.length,
        height=2*mavi_pad[2] + mavi_bbox.height,
    )

    start_end_walls = Configs.LinesCfg([
        [rect_cfg.bottom_left, rect_cfg.bottom_left + SVector(0, rect_cfg.height)],
        [rect_cfg.bottom_left + rect_cfg.size, rect_cfg.bottom_left + SVector(rect_cfg.length, 0)],
    ])

    space_cfg = Configs.SpaceCfg([
        (Configs.PeriodicWalls(), rect_cfg),
        (Configs.SlipperyWalls(), start_end_walls),
        (Configs.SlipperyWalls(), lines_cfg),
        [(Configs.SlipperyWalls(), c) for c in circles_list]...
    ])

    num_max_rings = 100
    rings_pos = Array{Float64, 3}(undef, 2, num_particles, num_max_rings)

    spawn_pos = create_circle((0, 0), ring_d/2, num_particles)
    spawn_pos = convert(Matrix{Float64}, spawn_pos)
    spawn_pos = copy(reinterpret(SVector{2, Float64}, vec(spawn_pos)))

    state = RingsState(
        rings_pos=rings_pos,
        pol=random_pol(num_max_rings),
        active_state=ActiveState(fill(false, num_max_rings)),
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
                bottom_left=rect_cfg.bottom_left + SVector(ring_d, rect_cfg.height/2 - ring_d),
                spawn_pos=spawn_pos,
                size=(2, 2),
                spawn_pol=:random,
                pad=get_particle_radius(dynamic_cfg),
            ),
            SinkCfg(Configs.RectangleCfg(
                length=ring_d*2,
                height=ring_d*2,
                bottom_left=rect_cfg.bottom_left + SVector(rect_cfg.length - 3*ring_d, rect_cfg.height/2 - ring_d),
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
            CircleGraphCfg(
                colors_map=:viridis,
            ),
            SinkSourceGraphCfg(),
        ])
    )
    
    # anim_cfg = VideoCfg(
    #     path="welcome_video.mp4",
    #     anim_cfg=anim_cfg,
    #     duration=20,
    # )

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