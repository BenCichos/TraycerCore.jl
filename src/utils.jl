function findroots(a::Float64, b::Float64, c::Float64)
    temp = b^2 - 4(a * c)
    temp < 1e-10 && return 0.0, 0.0
    return (-b - sqrt(temp)) / 2a, (-b + sqrt(temp)) / 2a
end

const AXIS2_X = @SVector [1, 0]
const AXIS2_Y = @SVector [0, 1]

const AXIS3_X = @SVector [1, 0, 0]
const AXIS3_Y = @SVector [0, 1, 0]
const AXIS3_Z = @SVector [0, 0, 1]

function getproperty(sv::SVector{3}, s::Symbol)
    s == :x && return sv[1]
    s == :y && return sv[2]
    s == :z && return sv[3]
    throw(ErrorException("type $(typeof(sv)) has no field $s"))
end

function getproperty(sv::SVector{2}, s::Symbol)
    s == :x && return sv[1]
    s == :y && return sv[2]
    throw(ErrorException("type $(typeof(sv)) has no field $s"))
end

Quaternion(v::SVector{2,<:Real}) = Quaternion(0, v[1], v[2], 0)
Quaternion(r::Real, v::SVector{2,<:Real}) = Quaternion(r, v[1], v[2], 0.0)
Quaternion(v::SVector{3,<:Real}) = Quaternion(0, v[1], v[2], v[3])
Quaternion(r::Real, v::SVector{3,<:Real}) = Quaternion(r, v[1], v[2], v[3])

SVector{2}(q::Quaternion) = SVector(q.v1, q.v2)
SVector{3}(q::Quaternion) = SVector(q.v1, q.v2, q.v3)

perpto(v::SVector{2,<:Real}) = SVector(v[2], -v[1])
perpto(v::SVector{3,<:Real}) = SVector(v[3], -v[2], v[1])

isparallel(u::SVector{N,<:Real}, v::SVector{N,<:Real}) where {N} = isone(dot(u, v))
isorthogonal(u::SVector{N,<:Real}, v::SVector{N,<:Real}) where {N} = iszero(dot(u, v))
isantiparallel(u::SVector{N,<:Real}, v::SVector{N,<:Real}) where {N} = isone(-dot(u, v))

*(q::Quaternion, v::SVector{N,<:Real}) where {N} = SVector{N}(q * Quaternion(v) * conj(q))

function quaternion(from::SVector{3,<:Real}, onto::SVector{3,<:Real})
    from, onto = normalize(from), normalize(onto)
    sin_theta, cos_theta = sincos(acos(dot(from, onto)) / 2)
    Quaternion(cos_theta, sin_theta * cross(from, onto))
end

function quaternion(axis::SVector{3,<:Real}, theta::Real, phi::Real)
    s_theta, c_theta = sincos(theta / 2)
    s_phi, c_phi = sincos(phi / 2)
    Quaternion(c_theta, s_theta * axis) * Quaternion(c_phi, s_phi * perpto(axis))
end

function quaternion(from::SVector{2,<:Real}, onto::SVector{2,<:Real})
    from, onto = normalize(from), normalize(onto)
    sin_theta, cos_theta = sincos(acos(dot(from, onto)) / 2)
    Quaternion(cos_theta, 0, 0, sin_theta)
end

rotate(vector::SVector{N,<:Real}, from::SVector{N,<:Real}, onto::SVector{N,<:Real}) where {N} = quaternion(from, onto) * vector
invrotate(vector::SVector{N,<:Real}, from::SVector{N,<:Real}, onto::SVector{N,<:Real}) where {N} = quaternion(onto, from) * vector
export rotate, invrotate

@inline quaternionx(from::SVector{2,<:Real}) = quaternion(from, AXIS2_X)
@inline quaterniony(from::SVector{2,<:Real}) = quaternion(from, AXIS2_Y)
@inline quaternionz(from::SVector{2,<:Real}) = quaternion(from, AXIS2_Y)

@inline quaternionx(from::SVector{3,<:Real}) = quaternion(from, AXIS3_X)
@inline quaterniony(from::SVector{3,<:Real}) = quaternion(from, AXIS3_Y)
@inline quaternionz(from::SVector{3,<:Real}) = quaternion(from, AXIS3_Z)

@inline invquaternionx(from::SVector{2,<:Real}) = quaternion(AXIS2_X, from)
@inline invquaterniony(from::SVector{2,<:Real}) = quaternion(AXIS2_Y, from)
@inline invquaternionz(from::SVector{2,<:Real}) = quaternion(AXIS2_Y, from)

@inline invquaternionx(from::SVector{3,<:Real}) = quaternion(AXIS3_X, from)
@inline invquaterniony(from::SVector{3,<:Real}) = quaternion(AXIS3_Y, from)
@inline invquaternionz(from::SVector{3,<:Real}) = quaternion(AXIS3_Z, from)

@inline rotatex(vector::SVector{N,<:Real}, from::SVector{N,<:Real}) where {N} = quaternionx(from) * vector
@inline rotatey(vector::SVector{N,<:Real}, from::SVector{N,<:Real}) where {N} = quaterniony(from) * vector
@inline rotatez(vector::SVector{N,<:Real}, from::SVector{N,<:Real}) where {N} = quaternionz(from) * vector

@inline invrotatex(vector::SVector{N,<:Real}, from::SVector{N,<:Real}) where {N} = invquaternionx(from) * vector
@inline invrotatey(vector::SVector{N,<:Real}, from::SVector{N,<:Real}) where {N} = invquaterniony(from) * vector
@inline invrotatez(vector::SVector{N,<:Real}, from::SVector{N,<:Real}) where {N} = invquaternionz(from) * vector

export findroots

export AXIS2_X, AXIS2_Y, AXIS3_X, AXIS3_Y, AXIS3_Z
export perpto, isparallel, isorthogonal, isantiparallel
export quaternion, quaternionx, quaterniony, quaternionz, invquaternionx, invquaterniony, invquaternionz
export rotate, invrotate, rotatex, rotatey, rotatez, invrotatex, invrotatey, invrotatez
