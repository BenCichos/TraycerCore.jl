module TraycerCore

using StaticArrays: SVector, SA, @SVector
using LinearAlgebra: normalize, norm, dot, cross
using Refraction: Material, isnullmaterial, NULL_MATERIAL, transmittance
import Refraction: transmittance
export material, isaligned

using Base: @kwdef
import Base: show, convert, push!, append!, iterate, pop!

include("ray.jl")
export AbstractRay
export Ray, TracedRay
export origin, direction, wavelength, index, external_indices, intensity, destination

include("generators.jl")
export SVectorGenerator, RayGenerator
export SVectors, Rays

include("interface.jl")
export AbstractInterface, Interface, CoatedInterface, ConsumingInterface
export isconsuminginterface, reflectance, transmittance

include("elements.jl")
export AbstractOpticalElement, PrimitiveElement, Object, Surface, CompoundElement, CompoundSurface, CompoundObject, CompoundObjectElement, CompoundSurfaceElement, CompoundElementPrimitive
export primitive_element, primitives, tracingtypes, normal, interface

include("opticalsystem.jl")
export OpticalSystem, rays, tracedrays, elements

include("tracingvector.jl")
export TracingVector, RecursiveTracingVector

include("intersection.jl")
export Intersection, intersection!, indices!, indices, incident, refraction, reflection, trigonometry, fresnel, snell
export point, reflect, refract

include("minintersection.jl")
export MinIntersection, distance!, minintersection!, update_index!, distance, intersection, index

include("trace.jl")
export trace

include("recursivetrace.jl")
export recursivetrace, _recursivetrace!


end
