abstract type AbstractOpticalElement{N} end

abstract type PrimitiveElement{N} <: AbstractOpticalElement{N} end
abstract type CompoundElement{N} <: AbstractOpticalElement{N} end

abstract type Object{N} <: PrimitiveElement{N} end
abstract type Surface{N} <: PrimitiveElement{N} end

abstract type CompoundObject{N} <: CompoundElement{N} end
abstract type CompoundSurface{N} <: CompoundElement{N} end

struct CompoundElementPrimitive{N,PE<:PrimitiveElement{N},O} <: PrimitiveElement{N}
    primitive_element::PE

    function CompoundElementPrimitive(primitive_element::PE, isobject::Bool) where {N,PE<:PrimitiveElement{N}}
        new{N,PE,isobject}(primitive_element)
    end
end
export CompoundElementPrimitive

const CompoundObjectElement{N,PE} = CompoundElementPrimitive{N,PE,true}
const CompoundSurfaceElement{N,PE} = CompoundElementPrimitive{N,PE,false}

primitive_element(compound_element::CompoundElementPrimitive) = compound_element.primitive_element

primitives(primitive_element::PE) where {PE<:PrimitiveElement} = [primitive_element,]
primitives(compoundelement::CE) where {CE<:CompoundElement} = [compoundelement.primitives...]
primitives(optical_elements::Vector{AOE}) where {N,AOE<:AbstractOpticalElement{N}} = Iterators.flatmap(primitives, optical_elements) |> collect

tracingtypes(primitive_element::PE) where {PE<:PrimitiveElement} = [primitive_element,]
tracingtypes(compound_object::CO) where {CO<:CompoundObject} = [CompoundElementPrimitive(primitive, true) for primitive in primitives(compound_object)]
tracingtypes(compound_surface::CS) where {CS<:CompoundSurface} = [CompoundElementPrimitive(primitive, false) for primitive in primitives(compound_surface)]
tracingtypes(optical_elements::Vector{<:AbstractOpticalElement{N}}) where {N} = Iterators.flatten([tracingtypes(optical_element) for optical_element in optical_elements])

interface(primitive_element::PE) where {PE<:PrimitiveElement} = primitive_element.interface
interface(compound_element_primitive::CompoundElementPrimitive) = interface(primitive_element(compound_element_primitive))
interface(element::AbstractOpticalElement{N}, ::SVector{N,Float64}) where {N} = interface(element)

isobject(object::Object) = true
isobject(compound_object::CompoundObject) = true
isobject(compound_object_element::CompoundObjectElement) = true

isobject(surface::Surface) = false
isobject(compound_surface::CompoundSurface) = false
isobject(compound_surface_element::CompoundSurfaceElement) = false

normal(element::AbstractOpticalElement{N}, ::SVector{N,Float64}) where {N} = normal(element)
normal(compound_element_primitive::CEP, point::SVector{N,Float64}) where {N,CEP<:CompoundElementPrimitive{N}} = normal(primitive_element(compound_element_primitive), point)
function normal(element::AbstractOpticalElement)
    throw(ErrorException("NotImplementedError: no method matching normal(:$(typeof(element))) has been implemented"))
    return nothing
end

onintersect(::AbstractOpticalElement{N}, ::Ray{N}, ::Float64, ::SVector{N,Float64}) where {N} = nothing
onintersect(compound_element_primitive::CEP, ray::Ray{N}, distance::Float64, normal::SVector{N,Float64}) where {N,CEP<:CompoundElementPrimitive{N}} = onintersect(primitive_element(compound_element_primitive), ray, distance, normal)
