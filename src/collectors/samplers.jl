abstract type SamplerCfg end
abstract type Sampler end

function get_sampler(cfg::SamplerCfg; kwargs...)
    throw(ArgumentError("not implemented for : $(typeof(cfg))"))
end

"""
    get_first_collect_index(cfg::SamplerCfg, ti)

Return the zero-based index of the first scheduled sample at or after `ti`.
"""
function get_first_collect_index(cfg::SamplerCfg, ti)
    throw(ArgumentError("not implemented for : $(typeof(cfg))"))
end

function get_num_collects(cfg::SamplerCfg, tf; ti=0)
    throw(ArgumentError("not implemented for : $(typeof(cfg))"))
end

function get_current_time(sampler::Sampler) 
    throw(ArgumentError("not implemented for : $(typeof(sampler))"))
end

update_sampler(sampler::Sampler, time_info) = update_sampler(sampler)
function update_sampler(sampler::Sampler) 
    throw(ArgumentError("not implemented for : $(typeof(sampler))"))
end

function to_collect(sampler::Sampler, time_info)
    throw(ArgumentError("not implemented for : $(typeof(sampler))"))
end

@kwdef struct LinearSamplerCfg{T} <: SamplerCfg
    to::T
    dt::T
    incremental::Bool=false
end

function get_sampler(cfg::LinearSamplerCfg; current_i=0)
    s = LinearSampler(cfg, current_i, zero(eltype(cfg.to)))
    s.next_time = get_current_time(s)
    return s
end

function get_first_collect_index(cfg::LinearSamplerCfg, ti)
    return max(0, ceil(Int, (ti - cfg.to) / cfg.dt))
end

function get_num_collects(cfg::LinearSamplerCfg, tf, ti=0)
    first_i = get_first_collect_index(cfg, ti)
    last_i = floor(Int, (tf - cfg.to) / cfg.dt)
    return max(0, last_i - first_i + 1)
end

mutable struct LinearSampler{T} <: Sampler
    cfg::LinearSamplerCfg{T}
    current_i::Int
    next_time::T
end
get_current_time(sampler::LinearSampler) = sampler.cfg.to + sampler.current_i * sampler.cfg.dt

function update_sampler(sampler::LinearSampler, time_info)
    if sampler.cfg.incremental
        sampler.next_time = time_info.time + sampler.cfg.dt
    else
        sampler.current_i += 1
        sampler.next_time = get_current_time(sampler)
    end
end

function to_collect(sampler::LinearSampler, time_info)
    return time_info.time >= sampler.next_time
end

# ===
# Log
# ===

@kwdef struct LogSamplerCfg{T} <: SamplerCfg
    xo::T
    dx::T
    base::T=10.0
    incremental::Bool=false
end

function get_sampler(cfg::LogSamplerCfg; current_i=0)
    s = LogSampler(cfg, current_i, zero(typeof(cfg.xo)))
    s.next_time = get_current_time(s)
    return s
end

function get_first_collect_index(cfg::LogSamplerCfg, ti)
    ti <= 0 && return 0
    return max(0, ceil(Int, (log(cfg.base, ti) - cfg.xo) / cfg.dx))
end

function get_num_collects(cfg::LogSamplerCfg, tf, ti=0)
    tf <= 0 && return 0

    first_i = get_first_collect_index(cfg, ti)
    last_i = floor(Int, (log(cfg.base, tf) - cfg.xo) / cfg.dx)

    return max(0, last_i - first_i + 1)
end

mutable struct LogSampler{T} <: Sampler
    cfg::LogSamplerCfg{T}
    current_i::Int
    next_time::T
end

get_current_time(sampler::LogSampler) = sampler.cfg.base^(sampler.cfg.xo + sampler.current_i * sampler.cfg.dx)

function update_sampler(sampler::LogSampler, time_info)
    if sampler.cfg.incremental
        sampler.next_time = time_info.time * sampler.cfg.base^sampler.cfg.dx
    else
        sampler.current_i += 1
        sampler.next_time = get_current_time(sampler)
    end
end

function to_collect(sampler::LogSampler, time_info)
    return time_info.time >= sampler.next_time
end
