
function trace(os::OpticalSystem{N}; depth::Int=10, threshold_intensity::Float64=1e-5) where {N}
    @assert depth >= 0 "Depth must be a non-negative integer"
    @assert threshold_intensity >= 0 "threshold_intensity must be a non-negative number"
    tracing_vector = TracingVector(rays(os), depth=depth, threshold_intensity=threshold_intensity)
    os_tracingtypes = AbstractOpticalElement{N}[tracingtypes(os)...]
    tracingtypes_indices = eachindex(os_tracingtypes)
    minintersection = MinIntersection()

    for current_depth in 1:depth
        for ray in tracing_vector
            reset!(minintersection)
            for index in tracingtypes_indices
                @inbounds minintersection!(minintersection, os_tracingtypes[index], ray)
                update_index!(minintersection, index)
            end
            isinf(distance(minintersection)) && continue
            @inbounds intersect!(tracing_vector, intersection(minintersection), os_tracingtypes[index(minintersection)], ray, distance(minintersection), current_depth == depth)
        end
        istraceablecounterzero(tracing_vector) && break
        step!(tracing_vector, current_depth)
    end
    traced_rays_array(tracing_vector)
end

trace!(os::OpticalSystem; depth::Int=10, threshold_intensity::Float64=1e-5) = (replace!(os, trace(os, depth=depth, threshold_intensity=threshold_intensity)); return os)
