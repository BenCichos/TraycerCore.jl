module TraycerCore

using Refraction: Material, isnullmaterial, NULL_MATERIAL, transmittance
export Material, transmittance
using StaticArrays: SVector, SA, @SVector
export SVector, SA, @SVector
using LinearAlgebra: normalize, norm, dot, cross
using Quaternions: Quaternion, imag_part, conj
using Base: @kwdef

using ApproximateRelations

@init_approx 1e-6

import Quaternions: Quaternion
import Refraction: transmittance
import Base: show, convert, push!, append!, iterate, :(*), getproperty


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

end
