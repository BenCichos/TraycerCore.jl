mutable struct Intersection{N}
    const indices::InterfaceIndices
    point::SVector{N,Float64}
    normal::SVector{N,Float64}
    isobject::Bool
end

Intersection() = Intersection(InterfaceIndices(), SVector{3,Float64}(NaN, NaN, NaN), SVector{3,Float64}(NaN, NaN, NaN), false)

function intersection!(intersection::Intersection, element::T, point::SVector{N,Float64}, isobject::Bool) where {N,T<:AbstractOpticalElement{N}}
    intersection.point = point
    intersection.normal = normal(element, point)
    intersection.isobject = isobject
end

point(intersection::Intersection) = intersection.point
normal(intersection::Intersection) = intersection.normal
indices(intersection::Intersection) = intersection.indices
incident(intersection::Intersection) = intersection.indices.incident
refraction(intersection::Intersection) = intersection.indices.refraction
reflection(intersection::Intersection) = intersection.indices.reflection
isobject(intersection::Intersection) = intersection.isobject

function intersect!(tracingvector::V, intersection::Intersection, element::T, ray::Ray{N}, distance::Float64, islastdepth::Bool) where {N,V<:AbstractTracingVector{N},T<:AbstractOpticalElement{N}}
    intersection!(intersection, element, origin(ray, distance), isobject(element))
    onintersect(element, ray, distance, normal(intersection))

    element_interface = interface(element, origin(ray, distance))
    islastdepth && return nothing
    isconsuminginterface(element_interface) && return nothing

    external_indices_vec = external_indices(ray) |> copy
    refraction_index = material(element_interface)(wavelength(ray))
    if isobject(intersection)
        if isaligned(ray, normal(intersection))
            refraction_index = pop!(external_indices_vec)
        else
            push!(external_indices_vec, refraction_index)
        end
    end

    indices!(indices(intersection), element_interface, index(ray), refraction_index, wavelength(ray))
    push!(tracingvector, TracedRay(ray, origin(ray, distance)))

    sin_incident, sin_transmitted, cos_incident, cos_transmitted = trigonometry(normal(intersection), direction(ray), incident(intersection), refraction(intersection))
    fresnel_reflectance, fresnel_transmittance = fresnel(incident(intersection), reflection(intersection), sin_incident, cos_incident, cos_transmitted)

    push!(tracingvector,
        Ray(
            origin(ray, distance),
            reflect(normal(intersection), direction(ray), cos_incident),
            intensity(ray) * fresnel_reflectance * reflectance(element_interface),
            wavelength(ray),
            incident(intersection),
            external_indices(ray),
        ))

    push!(tracingvector,
        Ray(
            origin(ray, distance),
            refract(normal(intersection), direction(ray), sin_incident, sin_transmitted, cos_incident, cos_transmitted),
            intensity(ray) * fresnel_transmittance * transmittance(element_interface, wavelength(ray), distance) * transmittance(element_interface),
            wavelength(ray),
            refraction(intersection),
            external_indices_vec,
        )
    )
    return nothing
end

function trigonometry(normal::SVector{N,Float64}, direction::SVector{N,Float64}, incident::Float64, transmitted::Float64) where {N}
    sin_incident, sin_transmitted = snell(normal, direction, incident, transmitted)
    cos_incident = dot(normal, direction)
    cos_transmitted = sqrt(1.0 - sin_transmitted^2)
    sin_incident, sin_transmitted, cos_incident, cos_transmitted
end

function snell(normal::SVector{N,Float64}, direction::SVector{N,Float64}, incident::Float64, transmitted::Float64) where {N}
    (normal == direction) || (normal == -direction) && return 0.0, 0.0
    sin_incident = norm(cross(normal, direction))
    sin_transmitted = rem(incident / transmitted * sin_incident, 1.0)
    sin_incident, sin_transmitted
end

function fresnel(incident::Float64, transmitted::Float64, sin_incident::Float64, cos_incident::Float64, cos_transmitted::Float64)
    sin_incident >= (transmitted / incident) && return 1.0, 0.0

    cos_incident = abs(cos_incident)
    cos_transmitted = abs(cos_transmitted)

    reflection_parallel = (incident * cos_transmitted - transmitted * cos_incident) / (incident * cos_transmitted + transmitted * cos_incident)
    reflection_perpendicular = (incident * cos_incident - transmitted * cos_transmitted) / (incident * cos_incident + transmitted * cos_transmitted)

    reflectance = 0.5 * (reflection_parallel^2 + reflection_perpendicular^2)
    transmittance = 1 - reflectance

    return reflectance, transmittance
end

function reflect(normal::SVector{N,Float64}, direction::SVector{N,Float64}, cos_incident::Float64) where {N}
    direction - 2 * normal * cos_incident
end

function refract(normal::SVector{N,Float64}, direction::SVector{N,Float64}, sin_incident::Float64, sin_transmitted::Float64, cos_incident::Float64, cos_transmitted::Float64) where {N}
    iszero(sin_incident) && return direction

    cos_incident = abs(cos_incident)
    cos_transmitted = abs(cos_transmitted)

    refracted_parallel = (direction - (cos_incident * normal)) * sin_transmitted / sin_incident
    refracted_perpendicular = normal * cos_transmitted
    normalize(refracted_parallel + refracted_perpendicular)
end
