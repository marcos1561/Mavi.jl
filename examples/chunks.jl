"""
Example: Using Chunks

This example creates many particles to demonstrate the difference
in performance between using chunks or not.
"""
module Example

using Mavi
using Mavi.Configs
using Mavi.InitStates
using Mavi.Visualization

using StaticArrays, StatsBase


function test()
    # Construct a SVector of Vectors
    # sv = SVector(Vector{Int}([1, 2, 3]), Vector{Int}([4, 5]), Vector{Int}())
    n = 1
    sv = SVector{n, Vector{SVector{2, Float64}}}([SVector{2, Float64}[(1, 2), (3, 4)] for _ in 1:n])
    
    for v in sv
        v .= Scalar(zero(eltype(v)))
    end

    a = [SVector(1.1, 2.2), SVector(1, 3)]
    # a[1] = a[1] .* SVector(2, -1)
    a[1]  *= SVector(2, -1)
    println(a)

end


function main()
    # Set to false to see the performance without chunks.
    use_chunks = true

    # Init state configs
    num_p_x = 25
    num_p_y = 25
    offset = 0.4

    num_p = num_p_x * num_p_y
    dynamic_cfg = HarmTruncCfg(10, 1, 1)
    radius = particle_radius(dynamic_cfg)

    pos, geometry_cfg = rectangular_grid(num_p_x, num_p_y, offset, radius)
    state = Mavi.SecondLawState(
        pos=pos,
        vel=random_vel(num_p, 1/5),
    )

    chunks_cfg = nothing
    if use_chunks    
        chunks_cfg = ChunksCfg(
            num_cols=trunc(Int, num_p_x*0.9), 
            num_rows=trunc(Int, num_p_y*0.9),
        )
    end
    
    int_cfg = IntCfg(
        dt=0.001,
        chunks_cfg=chunks_cfg,
    )

    system = System(
        state=state, 
        space_cfg=SpaceCfg(
            # wall_type=RigidWalls(),
            wall_type=PeriodicWalls(),
            geometry_cfg=geometry_cfg,
        ),
        dynamic_cfg=dynamic_cfg,
        int_cfg=int_cfg,
    )

    println(typeof(system.space_cfg))

    anim_cfg = Visualization.AnimationCfg(
        num_steps_per_frame=200,
        # num_steps_per_frame=50,
        exec_times_size=100,
        graph_cfg=CircleGraphCfg(colors_map=:magma),
    )

    # animate(system, Mavi.Integration.newton_step!, anim_cfg)
    
    function my_step!(system::System)
        clean_f = @timed Mavi.Integration.clean_forces!(system)
        chunk_up = @timed Mavi.Integration.update_chunks!(system.chunks)
        calc_f = @timed Mavi.Integration.calc_forces!(system)
        up = @timed Mavi.Integration.update_verlet!(system, Mavi.Integration.calc_forces!)
        walls = @timed Mavi.Integration.walls!(system, system.space_cfg)
        
        return clean_f.time, chunk_up.time, calc_f.time, up.time, walls.time
    end

    n = 200
    names = [
        :clean_f,
        :chunk_up,
        :calc_f,
        :up,
        :walls,
    ]

    
    buffers = (; (name => Vector{Float64}(undef, n) for name in names)...)

    my_step!(system)
    for i in 1:n
        times = my_step!(system)
        for (name, time) in zip(names, times)
            buffers[name][i] = time * 1000
        end
    end

    for name in names
        buf = buffers[name]
        m = round(mean(buf), digits=3)
        s = round(std(buf), digits=3)
        println("$name: $m ± $s")
    end

end

end

import .Example
Example.main()
# Example.test()