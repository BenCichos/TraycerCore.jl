abstract type SVectorGenerator{N} <: Function end

function SVectorGenerator(args::Union{AbstractRange,Real}...)
    ranges = filter(x -> x isa AbstractRange, args)
    @assert allequal(length.(ranges)) "All ranges must have the same length"
    N = length(first(ranges))
    converted_args = (isa(arg, AbstractRange) ? arg : range(arg, arg, N) for arg in args)
    Base.Generator(SVector, zip(converted_args...))
end

function SVectorGenerator{N}(args::Union{AbstractRange,Real}...) where {N}
    converted_args = map(args) do arg
        isa(arg, Real) && return range(arg, arg, N)
        range(arg[1], arg[end], length=N)
    end
    Base.Generator(SVector, zip(converted_args...))
end

function SVectorGenerator{N}(args::Vararg{Real}) where {N}
    converted_args = (range(arg, arg, length=N) for arg in args)
    Base.Generator(SVector, zip(converted_args...))
end

abstract type SVectors{N} <: Function end
SVectors(args::Union{AbstractRange,Real}...) = SVectorGenerator(args...) |> collect
SVectors{N}(args::Union{AbstractRange,Real}...) where {N} = SVectorGenerator{N}(args...) |> collect
SVectors{N}(args::Vararg{Real}) where {N} = SVectorGenerator{N}(args...) |> collect

function RayGenerator(origins::Base.Generator, directions::Base.Generator, intensity::Float64=1.0, wavelength::Float64=0.540)
    @assert length(origins) == length(directions) "The number of origins and directions must be the same"
    @assert length(origins.iter.is) == length(directions.iter.is) "The dimensions of origins and directions generator must be the same"
    N = length(origins.iter.is)
    (Ray(origins.f(vals[1:N]), directions.f(vals[N+1:end]), intensity, wavelength) for vals in zip(origins.iter.is..., directions.iter.is...))
end

function RayGenerator(origins::Vector{SVector{N,T}}, directions::Vector{SVector{N,T}}, intensity::Float64=1.0, wavelength::Float64=0.540) where {N,T<:Real}
    @assert !isempty(origins) "At least one ray must be provided"
    @assert length(origins) == length(directions) "The number of origins and directions must be the same"
    (Ray(origin, direction, intensity, wavelength) for (origin, direction) in zip(origins, directions))
end

function RayGenerator(origin::SVector{N,<:Real}, direction::SVector{N,<:Real}, intensities::AbstractRange, wavelengths::AbstractRange) where {N}
    @assert !isempty(intensities) "At least one ray must be provided"
    @assert length(intensities) == length(wavelengths) "The number of intensities and wavelengths must be the same"
    (Ray(origin, direction, intensity, wavelength) for (intensity, wavelength) in zip(intensities, wavelengths))
end

function RayGenerator(origins::Vector{SVector{N,T}}, directions::Vector{SVector{N,T}}, intensities::AbstractRange, wavelengths::AbstractRange) where {N,T<:Real}
    @assert !isempty(origins) "At least one ray must be provided"
    @assert length(origins) == length(directions) == length(intensities) == length(wavelengths) "The number of origins and directions must be the same"
    (Ray(origin, direction, intensity, wavelength) for (origin, direction, intensity, wavelength) in zip(origins, directions, intensities, wavelengths))
end

Rays(args...) = RayGenerator(args...) |> collect
