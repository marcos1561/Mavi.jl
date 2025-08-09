"Systems visualizations"
module Visualization

export animate, random_colors
export AnimationCfg, VideoCfg, UiSettings
export DefaultInfoUICfg
export MainGraphCfg, CircleGraphCfg, ScatterGraphCfg
export drawn_borders, colors_from_cmap

using GLMakie
using DataStructures

using Mavi.States: State
using Mavi.Systems
using Mavi.Configs
using Mavi.Utils.Progress

include("gui/info_ui.jl")
include("gui/system_graph.jl")

using .InfoUIs 
using .SystemGraphs 


function random_colors(num)
    return [rand(RGBf) for _ in 1:num]
end

"""
# Arguments
- sidebar_rel_length:  
    Relative horizontal space allocated to the InfoUI.
"""
@kwdef struct UiSettings
    sidebar_rel_length = 0.2
end

"""
Animation configurations

# Arguments
- graph_cfg:  
    Configuration for the part where the system is rendered.
    More info in "gui/system_graph.jl".

info_cfg:  
    Configuration for the information UI part.
    More info in "gui/info_ui.jl".

- fps:  
    Animation fps

- num_steps_per_frame:  
    How many time steps are done in a single frame.

- exec_times_size:  
    Circular buffer length that stores step time execution.

- ui_settings:  
    General settings for the UI window.
"""
@kwdef struct AnimationCfg{GraphT, InfoT}
    graph_cfg::GraphT = MainGraphCfg()
    info_cfg::InfoT = DefaultInfoUICfg()
    fps = 30
    num_steps_per_frame = 10
    exec_times_size = 100
    begin_paused = false
    ui_settings = UiSettings()
    fig_kwargs::Union{Dict, Nothing} = nothing
    ax_kwargs::Union{Dict, Nothing} = nothing
end

"""
Save video configuration.

- path: 
    path to save video (relative to where the REPL is)

- duration:
    Video duration in seconds.

- anim_cfg:
    Animation configurations.    
"""
struct VideoCfg{D<:Number, A<:AnimationCfg}
    path::String
    duration::D
    anim_cfg::A
end
function VideoCfg(;path, duration, anim_cfg=nothing)
    if anim_cfg === nothing
        anim_cfg = AnimationCfg()
    end

    VideoCfg(path, duration, anim_cfg)
end

"""
Information about the state of execution.

# Arguments
- sym_time:  
    Current simulation time (the physical time of the system).

- sym_time_count:  
    How many time steps were made.

- times:  
    Buffer with the execution time of the last performed time steps.   
"""
mutable struct ExecInfo 
    sym_time::Float64
    sym_time_count::Int
    times::CircularBuffer{Float64}
    times_ui::CircularBuffer{Float64}
end

get_anim_cfg(cfg::AnimationCfg) = cfg
get_anim_cfg(cfg::VideoCfg) = cfg.anim_cfg

get_graph_cfg(anim_cfg::AnimationCfg{G, I}) where {G<:GraphCfg, I} = anim_cfg.graph_cfg
get_graph_cfg(anim_cfg::AnimationCfg{G, I}) where {G<:SystemGraphs.GraphCompCfg, I} = MainGraphCfg(anim_cfg.graph_cfg)

"Render, in real time, the system using the given step function."
function animate(system::System, step!, cfg=nothing)
    GLMakie.activate!(; title="Mavi")

    if cfg === nothing
        cfg = AnimationCfg()
    end

    is_video = typeof(cfg) <: VideoCfg
    anim_cfg = get_anim_cfg(cfg)

    ax_kwargs = anim_cfg.ax_kwargs
    if ax_kwargs === nothing
        ax_kwargs = Dict()
    end

    fig_kwargs = anim_cfg.fig_kwargs
    if fig_kwargs === nothing
        fig_kwargs = Dict()
    end
    # if !haskey(fig_kwargs, :backgroundcolor)
    #     fig_kwargs[:backgroundcolor] = RGBf(0.94, 0.94, 0.94)
    # end

    fig = Figure(; fig_kwargs...)

    ui_settings = anim_cfg.ui_settings
    run_status = Observable(!anim_cfg.begin_paused)
    run_next_frame = false
    
    if !is_video
        system_gl = fig[1, 2] = GridLayout()
        system_ax = Axis(system_gl[1, 1]; aspect=DataAspect(), ax_kwargs...)
        # system_ax = Axis(system_gl[1, 1]; aspect=DataAspect())

        main_sidebar_gl = fig[1, 1] = GridLayout()

        sidebar_gl = main_sidebar_gl[1, 1] = GridLayout()
        Box(sidebar_gl[1, 1], cornerradius=5)
        Box(sidebar_gl[2, 1], cornerradius=5)
        
        info_gl = sidebar_gl[1, 1] = GridLayout(
            halign=:left, 
            valign=:top, 
            tellwidth=false,
        )
        
        control_gl = sidebar_gl[2, 1] = GridLayout(
            halign=:left, 
            valign=:top, 
            tellwidth=false,
        )
        rowsize!(control_gl, 1, Fixed(1))

        button_label = lift(run_status) do st
            if st
                return "Running"
            else 
                return "Stopped"
            end
        end
        start_stop_button = control_gl[2, 1] = Button(fig, 
            label=button_label,
            width=Relative(0.9),
        )
        on(start_stop_button.clicks) do n
            run_status[] = !run_status[]
            notify(run_status)
        end
        
        run_next_frame_button = control_gl[3, 1] = Button(fig, 
            label="Run Next Frame",
            width=Relative(0.9),
        )
        on(run_next_frame_button.clicks) do n
            if !run_status[]
                run_next_frame = true
            end   
        end

        # Box(main_sidebar_gl[:, 1], cornerradius=5)
        # Box(sidebar_gl[1, 1], color = (:blue, 0.1), strokecolor = :transparent)
        # Box(sidebar_gl[2, 1], color = (:red, 0.1), strokecolor = :transparent)
        Box(sidebar_gl[1, 1], cornerradius=5)
        # Box(sidebar_gl[2, 1], cornerradius=5)

        colsize!(fig.layout, 1, Relative(ui_settings.sidebar_rel_length))
        rowsize!(main_sidebar_gl, 1, Relative(1))
        
        rowsize!(sidebar_gl, 1, Relative(1/2))
        rowsize!(sidebar_gl, 2, Relative(1/2))
        
        
        info = InfoUIs.get_info_ui(info_gl, anim_cfg.info_cfg)
    else
        system_ax = Axis(fig[1, 1]; aspect=DataAspect(), ax_kwargs...)
    end

    graph = SystemGraphs.get_graph(system_ax, system, get_graph_cfg(anim_cfg))

    exec_info = ExecInfo(0, 0, 
        CircularBuffer{Float64}(anim_cfg.exec_times_size), 
        CircularBuffer{Float64}(anim_cfg.exec_times_size),
    )

    context = (
        anim_cfg=anim_cfg,
        system=system,
        graph=graph,
        exec_info=exec_info,
        run_status=true,
    )

    function make_frame(context)
        anim_cfg = context.anim_cfg
        system = context.system
        exec_info = context.exec_info
        graph = context.graph

        for _ in 1:anim_cfg.num_steps_per_frame
            step_info = @timed step!(system)
            push!(exec_info.times, step_info.time)
            exec_info.sym_time += system.int_cfg.dt
        end
        
        SystemGraphs.update_graph(graph, system)
    end
    
    if is_video
        num_frames = trunc(Int, anim_cfg.fps * cfg.duration)
        prog = ProgContinuos(init=1, final=num_frames)
        record(fig, cfg.path; framerate=anim_cfg.fps) do io
            for frame in 1:num_frames
                make_frame(context)
                recordframe!(io)
                show_progress(prog, frame)
            end
            show_finish(prog)
        end
    else
        display(fig)
        time_wait = 1/anim_cfg.fps 
        while isopen(fig.scene) 
            t1 = time()
            
            run_frame = false
            if run_status[]
                run_frame = true
            elseif run_next_frame
                run_next_frame = false
                run_frame = true
            end
            
            if run_frame
                make_frame(context)
            end
            
            InfoUIs.update_info_ui(info, exec_info, system)
            time_to_wait = time_wait - (time() - t1)
            if time_to_wait < 0
                time_to_wait = 0.01
            end
            sleep(time_wait/2)
            push!(exec_info.times_ui, time() - t1)
        end
    end
end

include("rings/view.jl")

end
