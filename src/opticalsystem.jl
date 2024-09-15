@kwdef struct OpticalSystem{N}
    elements::Vector{AbstractOpticalElement{N}}
    rays::Vector{Ray{N}} = Ray{N}[]
    tracedrays::Vector{TracedRay{N}} = TracedRay{N}[]

    function OpticalSystem(elements::Vector{AbstractOpticalElement{N}}, rays::Vector{Ray{N}}=Ray{N}[], tracedrays::Vector{TracedRay{N}}=TracedRay{N}[]) where {N}
        new{N}(elements, rays, tracedrays)
    end
end

elements(os::OpticalSystem) = os.elements
rays(os::OpticalSystem) = os.rays
tracedrays(os::OpticalSystem) = os.tracedrays
primitives(os::OpticalSystem) = primitives(elements(os))
tracingtypes(os::OpticalSystem) = tracingtypes(elements(os))

push!(os::OpticalSystem{N}, ray::Ray{N}) where {N} = push!(os.rays, ray)
push!(os::OpticalSystem{N}, traced_ray::TracedRay{N}) where {N} = push!(os.tracedrays, traced_ray)
push!(os::OpticalSystem{N}, element::AOE) where {N,AOE<:AbstractOpticalElement{N}} = push!(os.elements, element)

append!(os::OpticalSystem{N}, rays::Vector{Ray{N}}) where {N} = append!(os.rays, rays)
append!(os::OpticalSystem{N}, tracedrays::Vector{TracedRay{N}}) where {N} = append!(os.tracedrays, tracedrays)
append!(os::OpticalSystem{N}, elements::Vector{AOE}) where {N,AOE<:AbstractOpticalElement{N}} = append!(os.elements, elements)

reset!(os::OpticalSystem) = empty!(os.tracedrays)
empty!(os::OpticalSystem) = empty!(os.rays)
replace!(os::OpticalSystem{N}, rays::Vector{Ray{N}}) where {N} = (empty!(os); append!(os, rays))
replace!(os::OpticalSystem{N}, tracedrays::Vector{TracedRay{N}}) where {N} = (reset!(os); append!(os, tracedrays))

show(io::IO, os::OpticalSystem) = print(io,
    """
    OpticalSystem:
        $(length(elements(os))) elements
        $(length(rays(os))) rays
        $(length(tracedrays(os))) traced rays
    """
)


export OpticalSystem
export elements, rays, tracedrays, primitives, tracingtypes
export push!, append!, reset!, empty!, replace!
