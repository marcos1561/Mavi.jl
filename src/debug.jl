module Debug

using StaticArrays

export clean_debug!
export debug_pairwise_force!, debug_area_force!, debug_spring_force!
export DebugRingForces, ManyDebugs

abstract type DebugInfo end

function clean_debug!(debug_info) end

struct ManyDebugs{D}
    debugs::D
end
clean_debug!(debug_info::ManyDebugs) = clean_debug!.(debug_info.debugs)
    
end # Debug