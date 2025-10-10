module Progress
    using Printf

    export ProgContinuos, show_progress, show_finish
    export NormalFormatter, SilenceFormatter

    function seconds_to_hhmmss(seconds)
        h = trunc(Int, seconds ÷ 3600)
        m = trunc(Int, seconds ÷ 60 - 60 * h)
        s = trunc(Int, seconds) - 3600 * h - 60 * m 
        return lpad(h, 2, '0') * ":" * lpad(m, 2, '0') * ":" * lpad(s, 2, '0')
    end    

    abstract type ProgFormatter end

    struct SilenceFormatter <: ProgFormatter end

    function progress_text(formatter::SilenceFormatter, progress, remaining) end
    function final_text(prog, formatter::SilenceFormatter) end


    struct NormalFormatter{L} <: ProgFormatter
        label::L
    end

    function progress_text(formatter::NormalFormatter{Nothing}, progress, remaining)
        println("Progress: $(progress) % | $(remaining)")
    end
    function progress_text(formatter::NormalFormatter{String}, progress, remaining)
        l = formatter.label
        println("$l: Progress: $(progress) % | $(remaining)")
    end

    function final_text(prog, formatter::NormalFormatter{Nothing})
        elapsed_time = time() - prog.init_time
        println("Elapsed Time: $(seconds_to_hhmmss(elapsed_time))")
    end
    function final_text(prog, formatter::NormalFormatter{String})
        elapsed_time = time() - prog.init_time
        l = formatter.label
        println("$l: Elapsed Time: $(seconds_to_hhmmss(elapsed_time))")
    end

    mutable struct ProgInfo
        last_shown::Float64
        has_shown::Bool
    end

    struct ProgContinuos{T<:Number, F}
        init::T
        final::T
        show_step::Float64
        init_time::Float64
        check_show_time::Float64
        info::ProgInfo
        formatter::F
    end
    function ProgContinuos(;init, final, show_step=0.1, formatter=nothing)
        init, final = promote(init, final)

        if isnothing(formatter)
            formatter = NormalFormatter(nothing)
        end

        ProgContinuos(init, final, show_step, time(), 10.0, ProgInfo(0.0, false), formatter)
    end 

    function show_progress(prog, i)
        progress = (i - prog.init) / (prog.final - prog.init)
        
        info = prog.info
        to_show = (progress - info.last_shown) >= prog.show_step
        check_show = !info.has_shown && (time() - prog.init_time) > prog.check_show_time

        if to_show || check_show
            info.has_shown = true
            info.last_shown = progress

            elapsed = time() - prog.init_time
            
            if progress == 0
                remaining_txt = "∞"
            else
                remaining = elapsed * (1/progress - 1)
                remaining_txt = seconds_to_hhmmss(remaining)
            end
            
            progress_str = @sprintf("%05.2f", progress*100)

            progress_text(prog.formatter, progress_str, remaining_txt)
        end
    end

    show_finish(prog) = final_text(prog, prog.formatter)
end
