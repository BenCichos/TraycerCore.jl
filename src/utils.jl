Quaternion(v::SVector{3,<:Real}) = Quaternion(0, v[1], v[2], v[3])
Quaternion(r::Real, v::SVector{3,<:Real}) = Quaternion(r, v[1], v[2], v[3])

const _AXIS_X = @SVector [1, 0, 0]
const _AXIS_Y = @SVector [0, 1, 0]
const _AXIS_Z = @SVector [0, 0, 1]

perpendicular_vector(v::SVector{2,<:Real}) = SVector(v[2], -v[1])
perpendicular_vector(v::SVector{3,<:Real}) = SVector(v[3], -v[2], v[1])

isparallel(u::SVector{3,<:Real}, v::SVector{3,<:Real}) = isone(dot(u, v))
isorthogonal(u::SVector{3,<:Real}, v::SVector{3,<:Real}) = iszero(dot(u, v))
isantiparallel(u::SVector{3,<:Real}, v::SVector{3,<:Real}) = isone(-dot(u, v))

function rotation_from_zaxis(u::SVector{3,<:Real})
    u = normalize(u)
    norm = sqrt(u[1]^2 + u[2]^2)
    iszero(norm) && return (isone(u[3]) ? Quaternion(1, 0, 0, 0) : Quaternion(0, 1, 0, 0))
    s_theta, c_theta = sincos(acos(u[3]) / 2)
    s_phi, c_phi = sincos(acos(u[1] / norm) / 2)
    Quaternion(c_theta, u[2] * s_theta, -u[1] * s_theta, 0) * Quaternion(c_phi, s_phi * _AXIS_Z)
end

@inline rotation_to_zaxis(u::SVector{3,<:Real}) = inv(rotation_from_zaxis(u))

function align_normal_plane_with_zaxis(normal::SVector{3,<:Real}, vector::SVector{3,<:Real})
    inv(rotation_from_zaxis(normal)) * vector
end

function rotation_from_axis(axis::SVector{3,<:Real}, theta::Real)
    s, c = sincos(theta / 2)
    Quaternion(c, s * axis)
end

function rotation_from_axis(axis::SVector{3,<:Real}, theta::Real, phi::Real)
    s_theta, c_theta = sincos(theta / 2)
    s_phi, c_phi = sincos(phi / 2)
    Quaternion(c_theta, s_theta * axis) * Quaternion(c_phi, s_phi * perpendicular_vector(axis))
end

function rotate_vector(q::Quaternion, u::SVector{3,<:Real})
    SVector(imag_part(q * Quaternion(u) * conj(q)))
end
