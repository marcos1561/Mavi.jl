"""
Example: Using Rings with different types

Simulation with Rings using different ring types. Description:
- Two ring types: self-interactions have low adhesion, while pair interactions have high adhesion.
- The number of rings is constant.
- Rings are created in a regular rectangular grid with random polarizations and random types.
"""
module Example

using Mavi.Rings
using Mavi.Rings.Configs 
using Mavi.Visualization
using Mavi.Configs: get_bounding_box
using ColorSchemes

using Random
Random.seed!(1234)

const RGBf = ColorSchemes.RGB{Float32}

function create_system(;num_cols, num_rows)
    interaction_cfg_1 = HarmTruncCfg(
        k_rep = 13,
        k_atr = 0.1,
        dist_eq = 1,
        dist_max = 1 * 1.1,
    )

    interaction_cfg_2 = HarmTruncCfg(
        k_rep = 13,
        k_atr = 0.1,
        dist_eq = 0.8,
        dist_max = 0.8 * 1.1,
    )

    interaction_cfg_pair = HarmTruncCfg(
        k_rep = 13,
        k_atr = 3,
        dist_eq = 1,
        dist_max = 1 + 0.1,
    )
    
    interaction_finder = InteractionMatrix([
        [interaction_cfg_1, interaction_cfg_pair];; 
        [interaction_cfg_pair, interaction_cfg_2]
    ])

    interactions = list_interactions(interaction_finder)
    self_interactions = list_self_interactions(interaction_finder)

    dynamic_cfg = RingsCfg(
        p0=[3.5, 3.5],
        relax_time=[1., 1.],
        vo=[1., 1.],
        mobility=[1., 1.],
        rot_diff=[0.5, 0.5],
        k_area=[1., 1.],
        k_spring=[20., 20.],
        l_spring=[inter.dist_eq*0.8 for inter in self_interactions],
        num_particles=[10, 5],
        interaction_finder=interaction_finder,
    )

    interactions = list_interactions(dynamic_cfg.interaction_finder)

    num_rings = num_cols * num_rows
    types = rand(1:2, num_rings)
    p_radius = Configs.particle_radius.([interaction_cfg_1, interaction_cfg_2])
    rings_pos, geometry_cfg = InitStates.rectangular_grid(
        num_cols = num_cols,
        num_rows = num_rows,
        num_particles = dynamic_cfg.num_particles,
        p_radius = p_radius,
        types=types,
        pad_x=0.1, pad_y=0.1,
    )

    space_cfg = Configs.SpaceCfg(
        geometry_cfg = geometry_cfg,
        wall_type = Configs.PeriodicWalls(),
    )
    state = RingsState(
        rings_pos = rings_pos,
        pol = InitStates.random_pol(num_rings),
        types = types,
    )

    bbox = get_bounding_box(space_cfg.geometry_cfg)
    max_size = maximum([i.dist_max for i in interactions]) * 1.1
    num_x_chunks = floor(Int, bbox.length / max_size)
    num_y_chunks = floor(Int, bbox.height / max_size)

    system = RingsSystem(
        state = state,
        space_cfg = space_cfg,
        dynamic_cfg = dynamic_cfg,
        int_cfg = Configs.IntCfg(
            dt = 0.01,
            chunks_cfg = Configs.ChunksCfg(
                num_x_chunks, num_y_chunks,
            ),
        ),
    )

    return system
end

function main()
    system = create_system(
        num_cols = 6,
        num_rows = 6,
    )
    
    # colors_map = [:red, :blue]
    # colors::Vector{Symbol} = []
    # for ring_id in 1:size(system.state.rings_pos, 3)
    #     t = system.state.types[ring_id]
    #     push!(colors, colors_map[t])
    # end
    
    colors_map = ColorSchemes.bam
    # colors_map = ColorSchemes.vik
    colors::Vector{RGBf} = []
    for ring_id in 1:size(system.state.rings_pos, 3)
        t = system.state.types[ring_id]
        if t == 1
            c = rand(Float64)/4 + 3/4
        else    
            c = rand(Float64)/4
        end
        push!(colors, colors_map[c])
    end
    
    anim_cfg = AnimationCfg(
        num_steps_per_frame=10,
        graph_cfg=DefaultGraphCfg(colors_map=colors),
    )
    animate(system, Integration.step!, anim_cfg)

    # while true
    # for _ in 1:10
    #     Integration.step!(system, system.int_cfg)
    # end
end

end

import .Example
Example.main()