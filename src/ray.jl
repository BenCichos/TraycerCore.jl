abstract type AbstractRay{N} end

origin(ray::AR) where {AR<:AbstractRay} = ray.origin
intensity(ray::AR) where {AR<:AbstractRay} = ray.intensity
wavelength(ray::AR) where {AR<:AbstractRay} = ray.wavelength

@kwdef struct Ray{N} <: AbstractRay{N}
    origin::SVector{N,Float64}
    direction::SVector{N,Float64}
    intensity::Float64 = 1.0
    wavelength::Float64 = 0.54
    index::Float64 = 1.0
    external_indices::Vector{Float64} = Float64[]

    function Ray(origin::SVector{N,<:Real}, direction::SVector{N,<:Real}, intensity::Real=1.0, wavelength::Real=0.54, current_index::Real=1.0, external_indices::Vector{Float64}=Float64[]) where {N}
        @assert N == 2 || N == 3 "A ray must be either two- or three-dimensional"
        @assert 0 <= intensity <= 1 "Intensity must be between 0 and 1"
        @assert norm(direction) > 0 "Direction must be normalizable (non-zero)"
        new{N}(origin, normalize(direction), intensity, wavelength, current_index, external_indices)
    end


end

origin(ray::Ray, distance::Float64) = origin(ray) + distance * direction(ray)
direction(ray::Ray) = ray.direction
index(ray::Ray) = ray.index
external_indices(ray::Ray) = ray.external_indices
wavelength(ray::Ray) = ray.wavelength
wavelength(ray::Ray, newindex::Float64) = index(ray) * wavelength(ray) / newindex
propagate(ray::Ray, distance::Float64) = Ray(origin(ray) + distance * direction(ray), direction(ray), intensity(ray), wavelength(ray), index(ray), external_indices(ray))
isaligned(direction::SVector{N,Float64}, otherdirection::SVector{N,Float64}) where {N} = (dot(direction, otherdirection) > 0.0)
isaligned(ray::Ray{N}, otherdirection::SVector{N,Float64}) where {N} = isaligned(direction(ray), otherdirection)

show(io::IO, sys::Ray) = print(io, "Ray($(origin(sys)), $(direction(sys)), $(intensity(sys)), $(wavelength(sys)), $(index(sys)), $(external_indices(sys)))")

@kwdef struct TracedRay{N} <: AbstractRay{N}
    origin::SVector{N,Float64}
    destination::SVector{N,Float64}
    wavelength::Float64
    intensity::Float64
end

TracedRay(ray::Ray{N}, destination::SVector{N,Float64}) where {N} = TracedRay(origin(ray), destination, wavelength(ray), intensity(ray))
TracedRay(ray::Ray) = TracedRay(origin(ray), direction(ray), wavelength(ray), intensity(ray))

direction(traced_ray::TracedRay) = normalize(destination(traced_ray) - origin(traced_ray))
destination(traced_ray::TracedRay) = traced_ray.destination

show(io::IO, traced_ray::TracedRay) = print(io, "TracedRay($(origin(traced_ray)), $(destination(traced_ray)), $(wavelength(traced_ray)), $(intensity(traced_ray)))")
