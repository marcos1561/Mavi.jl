"Systems visualizations"
module Visualization

export animate, random_colors
export AnimationCfg, VideoCfg, ImageCfg, UiSettings
export DefaultInfoUICfg
export MainGraphCfg, CircleGraphCfg, ScatterGraphCfg, NumsGraphCfg
export drawn_borders, colors_from_cmap

using GLMakie
using DataStructures


using Mavi.States: State
using Mavi.Systems
using Mavi.Configs
using Mavi.MaviSerder
using Mavi.Utils.Progress
using Mavi.Integration: get_step_function

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
Animation configurations.

# Arguments
- graph_cfg:  
    Configuration for how the system is rendered.
    More info in "gui/system_graph.jl".

info_cfg:  
    Configuration for the information UI.
    More info in "gui/info_ui.jl".

- fps:  
    Animation fps

- num_steps_per_frame:  
    How many time steps are done in a single frame.

- slider_spf_range:
    Value range used by the slider to control `num_steps_per_frame`.

- exec_times_size:  
    Circular buffer length that stores step time execution.

- begin_paused:
    If `True` the animation starts paused.

- ui_settings:  
    General settings for the UI window.
"""
@kwdef struct AnimationCfg{GraphT, InfoT}
    graph_cfg::GraphT = MainGraphCfg()
    info_cfg::InfoT = DefaultInfoUICfg()
    fps = 30
    num_steps_per_frame = 10
    slider_spf_range = nothing
    exec_times_size = 100
    begin_paused = false
    ui_settings = UiSettings()
    fig_kwargs::Union{Dict, Nothing} = nothing
    ax_kwargs::Union{Dict, Nothing} = nothing
    save_fig_path::Union{String, Nothing} = "image.png"
end

"""
Save video configurations.

# Arguments
- path: 
    Path to save the video (relative to where the REPL is)

- duration:
    Video duration in seconds.

- anim_cfg:
    Animation configurations.    
"""
struct VideoCfg{D<:Number, A<:AnimationCfg}
    path::String
    duration::D
    anim_cfg::A
    save_configs::Bool
end
function VideoCfg(;path, duration, anim_cfg=nothing, save_configs=false)
    if anim_cfg === nothing
        anim_cfg = AnimationCfg()
    end

    valid_extensions = [".mp4", ".avi", ".mov", ".mkv", ".gif"]
    if !any(endswith(path, ext) for ext in valid_extensions)
        throw(ArgumentError("The path must point to a valid video file with one of the following extensions: $(join(valid_extensions, ", "))"))
    end

    mkpath(dirname(path))

    VideoCfg(path, duration, anim_cfg, save_configs)
end

@kwdef struct ImageCfg{GraphT, T<:Number}
    path::String
    tf::T = 0
    graph_cfg::GraphT = MainGraphCfg()
    fig_kwargs::Union{Dict, Nothing} = nothing
    ax_kwargs::Union{Dict, Nothing} = nothing
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
    times::CircularBuffer{Float64}
    times_ui::CircularBuffer{Float64}
end
ExecInfo(buffer_size) = ExecInfo(
    CircularBuffer{Float64}(buffer_size), 
    CircularBuffer{Float64}(buffer_size),
)

get_anim_cfg(cfg::AnimationCfg) = cfg
get_anim_cfg(cfg::VideoCfg) = cfg.anim_cfg

get_graph_cfg(anim_cfg::AnimationCfg{G, I}) where {G<:GraphCfg, I} = anim_cfg.graph_cfg
get_graph_cfg(anim_cfg::AnimationCfg{G, I}) where {G<:SystemGraphs.GraphCompCfg, I} = MainGraphCfg(anim_cfg.graph_cfg)
get_graph_cfg(image_cfg::ImageCfg{G, T}) where {G<:GraphCfg, T} = image_cfg.graph_cfg
get_graph_cfg(image_cfg::ImageCfg{G, T}) where {G<:SystemGraphs.GraphCompCfg, T} = MainGraphCfg(image_cfg.graph_cfg)

"Render, in real time, the system using the given step function."
function animate(system::System, cfg=nothing; step_func=nothing)
    if isnothing(step_func)
        step_func = get_step_function(system)
    end

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
    num_steps_per_frame = anim_cfg.num_steps_per_frame
    
    exec_info = ExecInfo(anim_cfg.exec_times_size)

    if !is_video
        system_gl = fig[1, 2] = GridLayout()
        system_ax = Axis(system_gl[1, 1]; aspect=DataAspect(), ax_kwargs...)
        # system_ax = Axis(system_gl[1, 1]; aspect=DataAspect())

        main_sidebar_gl = fig[1, 1] = GridLayout()

        sidebar_gl = main_sidebar_gl[1, 1] = GridLayout()
        
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

        Box(sidebar_gl[1, 1], cornerradius=5, color="#8d8d8dff")
        Box(sidebar_gl[2, 1], cornerradius=5, color="#8d8d8dff")

        button_label = lift(run_status) do st
            if st
                return "Running"
            else 
                return "Stopped"
            end
        end
        
        function start_stop_func()
            run_status[] = !run_status[]
            notify(run_status)
        end
        
        function run_next_frame_func()
            if !run_status[]
                run_next_frame = true
            end
        end

        function save_fig()
            path = anim_cfg.save_fig_path
            save(path, fig)
            println("Imagem Salva em: ", abspath(path))
        end

        start_stop_button = control_gl[2, 1] = Button(fig, 
            label=button_label,
            width=Relative(0.9),
        )
        on(start_stop_button.clicks) do n
            start_stop_func()
        end
        
        run_next_frame_button = control_gl[3, 1] = Button(fig, 
            label="Run Next Frame",
            width=Relative(0.9),
        )
        on(run_next_frame_button.clicks) do n
            run_next_frame_func()   
        end

        spf_range = anim_cfg.slider_spf_range
        num_steps_per_frame = anim_cfg.num_steps_per_frame
        if isnothing(spf_range)
            spf_range = 1:2*num_steps_per_frame
        end

        sliders = SliderGrid(
            control_gl[4, 1],
            (label = "Speed", range = spf_range, startvalue = num_steps_per_frame),
            # (label = "Current", range = 0:0.1:20, format = "{:.1f}A", startvalue = 10.2),
            # (label = "Resistance", range = 0:0.1:30, format = "{:.1f}Ω", startvalue = 15.9),
            width=Relative(0.9),
            # tellheight = true,
        )

        on(sliders.sliders[1].value) do speed
            num_steps_per_frame = speed
        end


        on(events(fig).keyboardbutton) do event
            if event.action == Keyboard.press
                key = event.key
                if key == Keyboard.right
                    run_next_frame_func()
                elseif key == Keyboard.space
                    start_stop_func()
                elseif key == Keyboard.h
                    reset_limits!(system_ax)
                elseif key == Keyboard.escape
                    GLMakie.close(display(fig))
                elseif key == Keyboard.s && 
                    Keyboard.left_control in events(fig).keyboardstate
                    save_fig()
                end
            end

             if event.action == Keyboard.press || event.action == Keyboard.repeat
                if event.key in (Keyboard.up, Keyboard.down)
                    step_size = 1
                    if Keyboard.left_shift in events(fig).keyboardstate
                        step_size = (maximum(spf_range) - minimum(spf_range))/5 
                    end
                    
                    step_sing = event.key == Keyboard.up ? 1 : -1 

                    slider = sliders.sliders[1]
                    set_close_to!(slider, num_steps_per_frame + step_sing*step_size)
                end
            end
        end

        # Box(main_sidebar_gl[:, 1], cornerradius=5)
        # Box(sidebar_gl[1, 1], color = (:blue, 0.1), strokecolor = :transparent)
        # Box(sidebar_gl[2, 1], color = (:red, 0.1), strokecolor = :transparent)
        # Box(sidebar_gl[1, 1], cornerradius=5)
        # Box(sidebar_gl[2, 1], cornerradius=5)

        colsize!(fig.layout, 1, Relative(ui_settings.sidebar_rel_length))
        rowsize!(main_sidebar_gl, 1, Relative(1))
        
        rowsize!(sidebar_gl, 1, Relative(1/2))
        rowsize!(sidebar_gl, 2, Relative(1/2))
        
        
        info = InfoUIs.get_info_ui(info_gl, anim_cfg.info_cfg)
        InfoUIs.update_info_ui(info, exec_info, system)
    else
        system_ax = Axis(fig[1, 1]; aspect=DataAspect(), ax_kwargs...)
    end

    graph = SystemGraphs.get_graph(system_ax, system, get_graph_cfg(anim_cfg))


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

        for _ in 1:num_steps_per_frame
            step_info = @timed step_func(system)
            push!(exec_info.times, step_info.time)
        end
        
        SystemGraphs.update_graph(graph, system)
    end

    function video_frame!(ax, system::System)
        time = round(system.time_info.time, digits=2)
        ax.title = "t = $time"
    end
    
    if is_video
        if cfg.save_configs
            parent_dir = dirname(cfg.path)
            filename = splitext(basename(cfg.path))[1]
            save_system_configs(system, parent_dir, "configs_$filename")
        end

        num_frames = trunc(Int, anim_cfg.fps * cfg.duration)
        prog = ProgContinuos(init=1, final=num_frames)
        record(fig, cfg.path; framerate=anim_cfg.fps) do io
            video_frame!(system_ax, system)
            recordframe!(io)
            for frame in 1:num_frames
                make_frame(context)
                video_frame!(system_ax, system)
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

function animate(system::System, cfg::ImageCfg, step! = nothing)
    if isnothing(step!)
        step! = get_step_function(system)
    end

    GLMakie.activate!(; title="Mavi")

    ax_kwargs = cfg.ax_kwargs
    if ax_kwargs === nothing
        ax_kwargs = Dict()
    end

    fig_kwargs = cfg.fig_kwargs
    if fig_kwargs === nothing
        fig_kwargs = Dict()
    end

    fig = Figure(; fig_kwargs...)
    ax = Axis(fig[1, 1]; aspect=DataAspect(), ax_kwargs...)

    graph = SystemGraphs.get_graph(ax, system, get_graph_cfg(cfg))

    ti = system.time_info.time
    while system.time_info.time - ti < cfg.tf
        step!(system)
    end

    SystemGraphs.update_graph(graph, system)

    save(cfg.path, fig)
end

# include("rings/view.jl")

end
