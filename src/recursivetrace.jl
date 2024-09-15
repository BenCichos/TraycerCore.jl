function _recursivetrace!(recursive_tracing_vector::RecursiveTracingVector{N}, minintersection::MinIntersection, ray::Ray{N}, depth::Int, threshold_intensity::Float64) where {N}
    (intensity(ray) <= threshold_intensity) && return nothing

    reset!(minintersection)
    for index in tracingtypes_indices(recursive_tracing_vector)
        @inbounds minintersection!(minintersection, tracingtypes(recursive_tracing_vector)[index], ray)
        update_index!(minintersection, index)
    end
    isinf(distance(minintersection)) && return nothing

    @inbounds intersect!(recursive_tracing_vector, intersection(minintersection), tracingtypes(recursive_tracing_vector)[index(minintersection)], ray, distance(minintersection), depth == 0)

    hasnewrays(recursive_tracing_vector) || return nothing

    reflected_ray, refracted_ray = rays(recursive_tracing_vector)
    _recursivetrace!(recursive_tracing_vector, minintersection, reflected_ray, depth - 1, threshold_intensity)
    _recursivetrace!(recursive_tracing_vector, minintersection, refracted_ray, depth - 1, threshold_intensity)

    return nothing
end

function recursivetrace(os::OpticalSystem{N}; depth::Int=10, threshold_intensity::Float64=1e-5, keeptracedrays::Bool=true) where {N}
    @assert depth >= 0 "depth must be a non-negative integer"
    @assert threshold_intensity >= 0 "threshold_intensity must be a non-negative number"

    recursive_tracing_vector = RecursiveTracingVector(2^3 * length(rays(os)), AbstractOpticalElement{N}[tracingtypes(os)...], keeptracedrays=keeptracedrays)
    minintersection = MinIntersection()

    for ray in rays(os)
        _recursivetrace!(recursive_tracing_vector, minintersection, ray, depth, threshold_intensity)
    end

    tracedrays(recursive_tracing_vector)
end
recursivetrace!(sys::OpticalSystem; depth::Int=10, threshold_intensity::Float64=1e-5) = (replace!(sys, recursivetrace(sys, depth=depth, threshold_intensity=threshold_intensity)); return sys)

export recursivetrace!, recursivetrace
