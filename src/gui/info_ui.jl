module InfoUIs

export InfoUICfg, DefaultInfoUICfg

using Printf
using GLMakie

abstract type InfoUICfg end

"""
Inform physical time (t) and average execution time per time step (Δt).
It's also possible to inject custom information as explained below.

# Arguments
- custom_items:  
    A function with the given signature:
        
        (System, ExecInfo) -> Array{Tuple{String, String}}
    
    Each tuple inside the returned array (item) should have the following information:
        item[1]: Item name. 
        item[2]: Item value (already formatted for printing).
        
    Items will be displayed line by line in the information panel in the following manner:
        item[1]: item[2]
"""
@kwdef struct DefaultInfoUICfg <: InfoUICfg 
    custom_items::Union{Function, Nothing} = nothing
end

struct DefaultInfoUI
    sym_time_obs::Observable{String}
    exec_times_obs::Observable{String}
    exec_times_ui_obs::Observable{String}
    custom_items::Observable{String}
    configs::DefaultInfoUICfg
end

"Build the InfoUI respective to `cfg` using the given `grid_layout` from Makie."
function get_info_ui(grid_layout, cfg::DefaultInfoUICfg)
    left_pad = 10
    up_pad = 10

    first_padding = (left_pad, 0, 0, up_pad)
    padding = (left_pad, 0, 0, 0)

    sym_time = Observable("t: 0")
    Label(grid_layout[1, 1], sym_time, halign=:left, padding=first_padding)

    exec_time = Observable("Δt: 0 ms")
    Label(grid_layout[2, 1], exec_time, halign=:left, padding=padding)
    
    exec_time_ui = Observable("FPS:")
    Label(grid_layout[3, 1], exec_time_ui, halign=:left, padding=padding)
    
    custom_items = Observable(" ")
    Label(grid_layout[4, 1], custom_items, halign=:left, padding=padding,
        justification=:left)

    DefaultInfoUI(sym_time, exec_time, exec_time_ui, custom_items, cfg)
end

function update_info_ui(info::DefaultInfoUI, exec_info, system)
    exec_times = exec_info.times
    exec_times_ui = exec_info.times_ui
    sym_time = system.time_info.time

    mean_exec_time = sum(exec_times)/length(exec_times) * 1000
    mean_exec_time_ui = sum(exec_times_ui)/length(exec_times_ui)
    info.sym_time_obs[] = @sprintf("t: %.3f", sym_time)
    info.exec_times_obs[] = @sprintf("Δt: %.3f ms", mean_exec_time)
    info.exec_times_ui_obs[] = @sprintf("FPS: %.2f", 1/mean_exec_time_ui)

    if !isnothing(info.configs.custom_items)
        items = info.configs.custom_items(system, exec_info)
    
        custom_string = ""
        for item in items
            custom_string *= "$(item[1]): $(item[2])\n"
        end
        info.custom_items[] = custom_string
    end
end

end