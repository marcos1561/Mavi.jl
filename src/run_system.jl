module RunSystem

export run_system

using Mavi.Integration: get_step_function, system_initialization
using Mavi.Utils.Progress

function run_system(system; tf=nothing, num_steps=nothing, step_func=nothing, init_func=nothing, stop_func=nothing)
    if isnothing(step_func)
        step_func = get_step_function(system)
    end
    
    if isnothing(init_func)
        system_initialization(system)
    else
        init_func(system)
    end

    if stop_func === nothing
        stop_func = _ -> false 
    end

    prog = nothing
    ti = system.time_info.time
    if !isnothing(tf)
        prog = ProgContinuos(init=system.time_info.time, final=tf)
        while system.time_info.time < tf
            step_func(system)
            show_progress(prog, system.time_info.time)
            if stop_func(system)
                break
            end
        end
    elseif !isnothing(num_steps)
        prog = ProgContinuos(init=system.time_info.num_steps, final=num_steps)
        while system.time_info.num_steps < num_steps
            step_func(system)
            show_progress(prog, count)
            if stop_func(system)
                break
            end 
        end
    else
        while true
            step_func(system)
            if stop_func(system)
                break
            end
        end
    end

    if prog !== nothing
        show_finish(prog)
    end
end

end # RunSystem