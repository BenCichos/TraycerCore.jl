abstract type AbstractInterface end

material(interface::I) where {I<:AbstractInterface} = interface.material
reflectance(interface::I) where {I<:AbstractInterface} = interface.reflectance
transmittance(interface::I) where {I<:AbstractInterface} = interface.transmittance
transmittance(interface::AbstractInterface, wavelength::Float64, distance::Float64) = transmittance(material(interface), wavelength, distance) * transmittance(interface)

@kwdef struct Interface <: AbstractInterface
    material::Material
    reflectance::Float64 = 1.0
    transmittance::Float64 = 1.0

    function Interface(material::Material, reflectance::Real=1.0, transmittance::Real=1.0)
        new(material, reflectance, transmittance)
    end
end

Interface(refractive_index::Real, reflectance::Real=1.0, transmittance::Real=1.0) = Interface(Material(refractive_index), reflectance, transmittance)

const NULL_INTERFACE = Interface(NULL_MATERIAL, 0.0, 1.0)
isnullinterface(interface::AbstractInterface) = (isnullmaterial(material(interface)) && iszero(reflectance(interface)) && isone(transmittance(interface)))

ConsumingInterface() = Interface(NULL_MATERIAL, 0.0, 0.0)
ConsumingInterface(m::Material) = Interface(m, 0.0, 0.0)
isconsuminginterface(interface::AbstractInterface) = iszero(reflectance(interface)) && iszero(transmittance(interface))

ReflectiveInterface(m::Material) = Interface(m, 1.0, 0.0)
isreflectiveinterface(interface::AbstractInterface) = isone(reflectance(interface)) && iszero(transmittance(interface))

TransmissiveInterface(m::Material) = Interface(m, 0.0, 1.0)
istransmissiveinterface(interface::AbstractInterface) = iszero(reflectance(interface)) && isone(transmittance(interface))

@kwdef struct CoatedInterface <: AbstractInterface
    elementinterface::Interface
    coatinginterface::Interface
    thickness::Float64
    reflectance::Float64 = 1.0
    transmittance::Float64 = 1.0

    function CoatedInterface(elementinterface::Interface, coatinginterface::Interface, thickness::Float64, reflectance::Real=1.0, transmittance::Real=1.0)
        new(elementinterface, coatinginterface, thickness, reflectance, transmittance)
    end
end

material(coatedinterface::CoatedInterface) = material(coatedinterface.elementinterface)
coatingmaterial(coatedinterface::CoatedInterface) = material(coatedinterface.coatinginterface)
thickness(coatedinterface::CoatedInterface) = coatedinterface.thickness
transmittance(coatedinterface::CoatedInterface, wavelength::Float64, ::Float64) = transmittance(coatingmaterial(coatedinterface), wavelength, thickness(coatedinterface)) * transmittance(coatedinterface)


convert(::Type{Interface}, m::Material) = Interface(m)
convert(::Type{Interface}, n::Real) = Interface(n)

mutable struct InterfaceIndices
    incident::Float64
    refraction::Float64
    reflection::Float64
end

InterfaceIndices() = InterfaceIndices(NaN, NaN, NaN)

function indices!(indices::InterfaceIndices, ::Interface, incident::Float64, refraction::Float64, ::Float64)
    indices.incident = incident
    indices.refraction = refraction
    indices.reflection = refraction
    return nothing
end

function indices!(indices::InterfaceIndices, interface::CoatedInterface, incident::Float64, refraction::Float64, wavelength::Float64)
    indices.incident = incident
    indices.refraction = refraction
    indices.reflection = coatingmaterial(interface)(wavelength)
    return nothing
end
