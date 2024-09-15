abstract type AbstractTracingVector{N} end

mutable struct TracingVector{N} <: AbstractTracingVector{N}
    const traceable_rays_data::Vector{Ray{N}}
    const traced_rays_data::Vector{TracedRay{N}}
    const threshold_intensity::Float64
    const depth::Int
    current_allocated_depth::Int
    previous_traceable_length::Int
    current_traceable_range::UnitRange{Int}
    traceable_rays_counter::Int
    traced_rays_counter::Int
end

function TracingVector(rays::Vector{Ray{N}}; depth::Int, threshold_intensity::Float64) where {N}
    current_allocated_depth = depth <= 3 ? depth : 3
    traceable_rays_data = Vector{Ray{N}}(undef, 2^current_allocated_depth * length(rays))
    traced_rays_data = Vector{TracedRay{N}}(undef, 2^current_allocated_depth * length(rays))

    @inbounds for (i, ray) in enumerate(rays)
        traceable_rays_data[i] = ray
    end

    TracingVector{N}(
        traceable_rays_data,
        traced_rays_data,
        threshold_intensity,
        depth,
        current_allocated_depth,
        length(traceable_rays_data),
        1:length(rays),
        0,
        0
    )
end

traceable_rays_data(v::TracingVector) = v.traceable_rays_data
traced_rays_data(v::TracingVector) = v.traced_rays_data
threshold_intensity(v::TracingVector) = v.threshold_intensity
traceable_rays_counter(v::TracingVector) = v.traceable_rays_counter
traced_rays_counter(v::TracingVector) = v.traced_rays_counter
length_traceable_rays(v::TracingVector) = length(v.traceable_rays_data)
length_traced_rays(v::TracingVector) = length(v.traced_rays_data)
istraceablecounterzero(v::TracingVector) = iszero(traceable_rays_counter(v))
current_traceable_range(v::TracingVector) = v.current_traceable_range
depth(v::TracingVector) = v.depth
current_allocated_depth(v::TracingVector) = v.current_allocated_depth
previous_traceable_length(v::TracingVector) = v.previous_traceable_length
traced_rays_array(v::TracingVector) = v.traced_rays_data[1:traced_rays_counter(v)]

increment_traceable!(v::TracingVector) = v.traceable_rays_counter += 1
increment_traced!(v::TracingVector) = v.traced_rays_counter += 1
reset_traceable_counter!(v::TracingVector) = v.traceable_rays_counter = 0
current_traceable_range!(v::TracingVector, new_range::UnitRange{Int}) = (v.current_traceable_range = new_range)
function extend!(v::TracingVector, n::Int)
    resize!(traceable_rays_data(v), traceable_rays_counter(v) + n)
    resize!(traced_rays_data(v), traced_rays_counter(v) + n)
    v
end
current_allocated_depth!(v::TracingVector, new_allocated_depth::Int) = (v.current_allocated_depth = new_allocated_depth)
previous_traceable_length!(v::TracingVector) = (v.previous_traceable_length = length_traceable_rays(v))

iterate(v::TracingVector, state::Int=1) = state > length(current_traceable_range(v)) ? nothing : (getindex(v, current_traceable_range(v)[state]), state + 1)
getindex(v::TracingVector, i::Int) = @inbounds traceable_rays_data(v)[mod1(i, previous_traceable_length(v))]

function push!(v::TracingVector{N}, ray::Ray{N}) where {N}
    intensity(ray) <= threshold_intensity(v) && return nothing
    @inbounds traceable_rays_data(v)[mod1(increment_traceable!(v) + last(current_traceable_range(v)), length_traceable_rays(v))] = ray
    nothing
end

function push!(v::TracingVector{N}, reflected_ray::Ray{N}, refracted_ray::Ray{N}) where {N}
    if intensity(reflected_ray) <= threshold_intensity(v)
        @inbounds traceable_rays_data(v)[mod1(increment_traceable!(v) + last(current_traceable_range(v)), length_traceable_rays(v))] = reflected_ray
    end
    if intensity(refracted_ray) <= threshold_intensity(v)
        @inbounds traceable_rays_data(v)[mod1(increment_traceable!(v) + last(current_traceable_range(v)), length_traceable_rays(v))] = refracted_ray
    end
    nothing
end

function push!(v::TracingVector{N}, traced_ray::TracedRay{N}) where {N}
    @inbounds traced_rays_data(v)[increment_traced!(v)] = traced_ray
    v
end

function step!(tracing_vector::TracingVector, current_depth::Int)
    previous_traceable_length!(tracing_vector)
    old_range = current_traceable_range(tracing_vector)
    start_index = last(old_range) + 1
    last_index = last(old_range) + traceable_rays_counter(tracing_vector)
    new_range = start_index:last_index
    if current_depth == current_allocated_depth(tracing_vector)
        new_allocated_depth = current_allocated_depth(tracing_vector) + 3
        adjustment = new_allocated_depth <= depth(tracing_vector) ? 0 : (depth(tracing_vector) - new_allocated_depth)
        extended_allocated_depth = 3 + adjustment
        current_allocated_depth!(tracing_vector, new_allocated_depth + adjustment)

        extend!(tracing_vector, 2^extended_allocated_depth * traceable_rays_counter(tracing_vector))
    end
    current_traceable_range!(tracing_vector, new_range)
    reset_traceable_counter!(tracing_vector)
end

mutable struct RecursiveTracingVector{N} <: AbstractTracingVector{N}
    const tracedrays::Vector{TracedRay{N}}
    const tracingtypes::Vector{AbstractOpticalElement{N}}
    const tracingtypes_indices::Base.OneTo{Int64}
    const traceablerays::Vector{Ray{N}}
    const keeptracedrays::Bool
    hasnewrays::Bool

    function RecursiveTracingVector(sizehint::Int, tracingtypes::Vector{AbstractOpticalElement{N}}; keeptracedrays::Bool=true) where {N}
        traced_rays_array = Vector{TracedRay{N}}()
        keeptracedrays && sizehint!(traced_rays_array, sizehint)
        traceablerays = Vector{Ray{N}}(undef, 2)
        new{N}(traced_rays_array, tracingtypes, eachindex(tracingtypes), traceablerays, keeptracedrays, false)
    end
end

tracedrays(vector::RecursiveTracingVector) = vector.tracedrays
tracingtypes(vector::RecursiveTracingVector) = vector.tracingtypes
tracingtypes_indices(vector::RecursiveTracingVector) = vector.tracingtypes_indices
traceablerays(vector::RecursiveTracingVector) = vector.traceablerays
hasnewrays(vector::RecursiveTracingVector) = vector.hasnewrays
hasnewrays!(vector::RecursiveTracingVector, newrays::Bool) = (vector.hasnewrays = newrays)
keeptracedrays(vector::RecursiveTracingVector) = vector.keeptracedrays

function push!(vector::RecursiveTracingVector{N}, ray::Ray{N}) where {N}
    hasnewrays(vector) && (traceablerays(vector)[2] = ray; return nothing)
    hasnewrays!(vector, true)
    traceablerays(vector)[1] = ray
    nothing
end

function rays(vector::RecursiveTracingVector{N}) where {N}
    hasnewrays!(vector, false)
    return traceablerays(vector)[1], traceablerays(vector)[2]
end

function push!(vector::RecursiveTracingVector{N}, traced_ray::TracedRay{N}) where {N}
    keeptracedrays(vector) || return nothing
    push!(tracedrays(vector), traced_ray)
    return nothing
end
