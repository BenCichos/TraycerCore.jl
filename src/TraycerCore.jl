module TraycerCore

using Refraction: Material, isnullmaterial, NULL_MATERIAL, transmittance
using StaticArrays: SVector, SA, @SVector
using LinearAlgebra: normalize, norm, dot, cross
using Quaternions: Quaternion, imag_part, conj
using Base: @kwdef
using ApproximateRelations

import Quaternions: Quaternion
import Refraction: transmittance
import Base: show, convert, push!, append!, iterate, :(*)

export SVector, SA, @SVector
export Material, NULL_MATERIAL, isnullmaterial, transmittance

include("utils.jl")

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
include("box.jl")

__init__() = set_approx!(1e-6)

end
