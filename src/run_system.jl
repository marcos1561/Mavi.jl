module RunSystem

export run_system

using Mavi.Integration: get_step_function
using Mavi.Utils.Progress

function run_system(system; tf=nothing, num_steps=nothing, step_func=nothing)
    if isnothing(step_func)
        step_func = get_step_function(system)
    end
    prog = nothing
    ti = system.time_info.time
    if !isnothing(tf)
        prog = ProgContinuos(init=system.time_info.time, final=tf)
        while system.time_info.time < tf
            step_func(system)
            show_progress(prog, system.time_info.time)
        end
    else
        prog = ProgContinuos(init=system.time_info.num_steps, final=num_steps)
        while system.time_info.num_steps < num_steps
            step_func(system)
            show_progress(prog, count)
        end
    end
    show_finish(prog)
end

end # RunSystem