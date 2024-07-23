"Systems visualizaions"
module Visualization

using GLMakie, Printf
using Mavi: State, System
using Mavi.Configs
using DataStructures

include("gui/info_ui.jl")
include("gui/system_graph.jl")

using .InfoUI 
using .SystemGraphs 

@kwdef struct UiSettings
    sidebar_rel_length = 0.2
end

"""
Animation configurations

# Arguments
- fps: 
Animation fps
- num_stesp_per_frame: 
    How many time steps are done in a single frame.
- circle_radius: 
    Radius used to draw circles centered on particles positions.  
    If not given, `particle_radius(dynamic_cfg)` will be used.
- circle_rel: 
    How many verticies used to draw circles.
- exec_times_size: 
    Circular buffer length that stores step time execution.
"""
@kwdef struct AnimationCfg
    fps = 30
    num_steps_per_frame = 10
    circle_radius = -1.0
    circle_rel = 20
    exec_times_size = 40
    ui_settings = UiSettings()
end


"Return circle vertices centered in the origin with given radius."
function get_circle(radius, resolution=20)
    theta = LinRange(0, 2π, resolution+1)

    circle = Vector{Point2f}()
    for id in 1:length(theta)-1
        t1, t2 = theta[id], theta[id+1]
        push!(circle, Point2f(radius * cos(t1), radius * sin(t1))) 
        push!(circle, Point2f(radius * cos(t2), radius * sin(t2))) 
    end

    return circle
end

function drawn_borders(ax, space_cfg::RectangleCfg)
    l, h = space_cfg.length, space_cfg.height
    lines!(ax, [0, l], [0, 0], color=:black)
    lines!(ax, [0, l], [h, h], color=:black)
    lines!(ax, [0, 0], [0, h], color=:black)
    lines!(ax, [l, l], [0, h], color=:black)
end

function drawn_borders(ax, space_cfg::CircleCfg)
    arc!(ax, Point2f(0), space_cfg.radius, -π, π, color=:black)
end

function update_circles!(circle_points, circle_base, state::State{T}) where {T}
    num_edges = length(circle_base)
    num_p = size(state.pos)[2]
    for i in 1:num_p
        begin_id = num_edges * (i-1) + 1
        end_id = begin_id + num_edges - 1 
        center = [state.pos[1, i], state.pos[2, i]]
        circle_points[begin_id:end_id] .= circle_base .+ [center]
    end
end

mutable struct ExecInfo 
    sym_time::Float64
    times::CircularBuffer{Float64}
end

"Render, in real time, the system using the given step function."
function animate(system::System{T}, step!, cfg::AnimationCfg) where {T}
    x_obs = Observable(@view system.state.pos[1, :])
    y_obs = Observable(@view system.state.pos[2, :])
    pos_obs = [x_obs, y_obs]
    
    fig = Figure(
        backgroundcolor = RGBf(0.98, 0.98, 0.98), 
        # size = (1000, 700),
    )

    system_gl = fig[1, 2] = GridLayout()
    
    sidebar_gl = fig[1, 1] = GridLayout()
    info_gl = sidebar_gl[1, 1] = GridLayout(
        halign=:left, 
        valign=:top, 
        tellwidth=false,
    )
    
    rowsize!(sidebar_gl, 1, Relative(1))
    colsize!(fig.layout, 1, Relative(0.2))
    
    Box(sidebar_gl[:, 1], cornerradius=5)


    label_text = Observable("t: 0")
    Label(info_gl[1, 1], label_text)

    ax = Axis(system_gl[1, 1])
    scatter!(ax, pos_obs...)
    # ax.aspect = 1
    
    drawn_borders(ax, system.space_cfg)

    if cfg.circle_radius == -1
        radius = particle_radius(system.dynamic_cfg) 
    else
        radius = cfg.circle_radius 
    end
    
    # Circles
    circle_base = get_circle(radius, cfg.circle_rel)
    circle_points = Vector{Point2f}(undef, system.num_p * length(circle_base))
    
    update_circles!(circle_points, circle_base, system.state)
    circle_points_obs = Observable(circle_points)
    circles_lines = linesegments!(ax, circle_points_obs, color=:black)

    exec_info = ExecInfo(0, CircularBuffer{Float64}(cfg.exec_times_size))

    display(fig)
    while events(fig).window_open[] 
        for _ in 1:cfg.num_steps_per_frame
            info = @timed step!(system, system.int_cfg)
            push!(exec_info.times, info.time)
            exec_info.sym_time += system.int_cfg.dt
        end
        
        mean_exec_time = sum(exec_times)/length(exec_times) * 1000
        ax.title = @sprintf("t=%.3f | Δt=%.3f ms", time, mean_exec_time)

        label_text[] = @sprintf("t: %.3f", time)

        update_circles!(circle_points, circle_base, system.state)
        notify(circle_points_obs)

        notify(pos_obs[1])
        notify(pos_obs[2])
        sleep(1/cfg.fps)
    end
end

end
