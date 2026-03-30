module RunSystem

export run_system

using Mavi.Integration: get_step_function
using Mavi.Utils.Progress

function run_system(system; tf=nothing, num_steps=nothing, step_func=nothing)
    if isnothing(step_func)
        step_func = get_step_function(system)
    end
    prog = nothing
    if !isnothing(tf)
        prog = ProgContinuos(init=0, final=tf)
        while system.time_info.time < tf
            step_func(system)
            show_progress(prog, system.time_info.time)
        end
    else
        count = 0
        prog = ProgContinuos(init=0, final=num_steps)
        while count < num_steps
            step_func(system)
            count += 1
            show_progress(prog, count)
        end
    end
    show_finish(prog)
end

end # RunSystem