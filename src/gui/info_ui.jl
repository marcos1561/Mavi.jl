module InfoUIs

export InfoUICfg, DefaultInfoUICfg

using Printf
using GLMakie

abstract type InfoUICfg end

struct DefaultInfoUICfg <: InfoUICfg end

struct DefaultInfoUI
    sym_time_obs::Observable{String}
    exec_times_obs::Observable{String}
end

function get_info_ui(grid_layout, cfg::DefaultInfoUICfg)
    left_pad = 10
    up_pad = 10

    first_padding = (left_pad, 0, 0, up_pad)
    padding = (left_pad, 0, 0, 0)

    sym_time = Observable("t: 0")
    Label(grid_layout[1, 1], sym_time, padding=first_padding)
   
    exec_time = Observable("Δt: 0 ms")
    Label(grid_layout[2, 1], exec_time, padding=padding)
    
    DefaultInfoUI(sym_time, exec_time)
end

function update_info_ui(info::DefaultInfoUI, exec_info)
    exec_times = exec_info.times
    sym_time = exec_info.sym_time

    mean_exec_time = sum(exec_times)/length(exec_times) * 1000
    info.sym_time_obs[] = @sprintf("t: %.3f", sym_time)
    info.exec_times_obs[] = @sprintf("Δt: %.3f ms", mean_exec_time)
end

end