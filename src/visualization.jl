"Systems visualizaions"
module Visualization

using GLMakie, Printf
using Mavi: State, System
using Mavi.Configs

@kwdef struct AnimationCfg
    fps::Int=30
    num_steps_per_frame::Int
    circle_radius::Float64=-1
    circle_rel::Int = 20
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
    num_p = length(state.x)
    for i in 1:num_p
        begin_id = num_edges * (i-1) + 1
        end_id = begin_id + num_edges - 1 
        center = [state.x[i], state.y[i]]
        circle_points[begin_id:end_id] .= circle_base .+ [center]
    end
end

"Render, in real time, the system using the given step function."
function animate(system::System{T}, step!, cfg::AnimationCfg) where {T}
    x_obs = Observable(system.state.x)
    y_obs = Observable(system.state.y)
    pos_obs = [x_obs, y_obs]
    
    fig, ax = scatter(pos_obs...)
    ax.aspect = 1
    
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

    display(fig)
    time = 0.0
    while events(fig).window_open[] 
        for _ in 1:cfg.num_steps_per_frame
            step!(system)
            time += system.int_cfg.dt
        end
        
        ax.title = @sprintf("t=%.3f", time)

        update_circles!(circle_points, circle_base, system.state)
        notify(circle_points_obs)

        notify(pos_obs[1])
        notify(pos_obs[2])
        sleep(1/cfg.fps)
    end
end

end
