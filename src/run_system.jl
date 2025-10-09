module RunSystem

export run_system

using Mavi.Integration: get_step_function

function run_system(system; tf=nothing, num_steps=nothing, step_func=nothing)
    if isnothing(step_func)
        step_func = get_step_function(system)
    end

    if !isnothing(tf)
        while system.time_info.time < tf
            step_func(system)
        end
    else
        count = 0
        while count < num_steps
            step_func(system)
            count += 1
        end
    end
end

end # RunSystem