module TraycerCore

using StaticArrays: SVector, SA, @SVector
using LinearAlgebra: normalize, norm, dot, cross
using Refraction: Material, isnullmaterial, NULL_MATERIAL, transmittance
import Refraction: transmittance

using Base: @kwdef
import Base: show, convert, push!, append!, iterate

export SVector, SA, @SVector
export normalize, norm, dot, cross
export Material, NULL_MATERIAL, isnullmaterial, transmittance

include("ray.jl")
include("generators.jl")
include("interface.jl")
include("elements.jl")
include("opticalsystem.jl")
include("tracingvector.jl")
include("intersection.jl")
include("minintersection.jl")
include("trace.jl")
include("recursivetrace.jl")


end
