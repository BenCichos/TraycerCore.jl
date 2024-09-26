@kwdef struct Box{N} <: Object{N}
    origin::SVector{N,Float64}
    size::SVector{N,Float64}
    axis::SVector{N,Float64}
    interface::Interface

    function Box(origin::SVector{2,<:Real}, size::SVector{2,<:Real}, axis::SVector{2,<:Real}, interface=NULL_INTERFACE)
        new{2}(origin, size, axis, interface)
    end

    function Box(origin::SVector{3,<:Real}, size::SVector{3,<:Real}, axis::SVector{3,<:Real}, interface=NULL_INTERFACE)
        new{3}(origin, size, axis, interface)
    end
end
export Box

origin(box::Box) = box.origin
size(box::Box) = box.size
axis(box::Box) = box.axis
bounds(box::Box) = (origin(box), origin(box) + size(box))
hassize(box::Box) = !all(iszero, size(box))

@approx function isinside(box::Box{N}, vector::SVector{N,Float64}) where {N}
    hassize(box) || return false
    axis_vector = invquaterniony(axis(box)) * (vector - origin(box))

    for point in axis_vector
        point < 0.0 && return false
    end

    for (point, axis_size) in zip(axis_vector, size(box))
        point > axis_size && return false
    end

    return true
end

@approx function onsurface(box::Box{N}, vector::SVector{N,Float64}) where {N}
    axis_vector = invquaterniony(axis(box)) * (vector - origin(box))
    for point in axis_vector
        point < 0.0 && return false
    end

    for (point, axis_size) in zip(axis_vector, size(box))
        point > axis_size && return false
    end

    for (point, axis_size) in zip(axis_vector, size(box))
        point == 0 || point == axis_size && return true
    end

    return false
end

@approx function normal(box::Box{2}, vector::SVector{2,Float64})
    x, y = size(box)
    axis_vector = invquaterniony(axis(box)) * (vector - origin(box))
    normal_x, normal_y = 0.0, 0.0
    if axis_vector.x == 0.0
        normal_x = -1.0
    elseif axis_vector[1] == x
        normal_x = 1.0
    end
    if axis_vector[2] == 0.0
        normal_y = -1.0
    elseif axis_vector[2] == y
        normal_y = 1.0
    end

    quaterniony(axis(box)) * SA[normal_x, normal_y]
end

@approx function normal(box::Box{3}, vector::SVector{3,Float64})
    x, y, z = size(box)
    axis_vector = invquaterniony(axis(box)) * (vector - origin(box))
    normal_x, normal_y, normal_z = 0.0, 0.0, 0.0
    if axis_vector.x == 0.0
        normal_x = -1.0
    elseif axis_vector.x == x
        normal_x = 1.0
    end

    if axis_vector.y == 0.0
        normal_y = -1.0
    elseif axis_vector.y == y
        normal_y = 1.0
    end

    if axis_vector.z == 0.0
        normal_y = -1.0
    elseif axis_vector.z == z
        normal_y = 1.0
    end

    quaternionz(axis(box)) * SA[normal_x, normal_y, normal_z]
end

@approx function minintersection!(minintersection::MinIntersection, box::Box{3}, ray::Ray{3})
    hassize(box) || return nothing
    box_invquaterion = invquaternionz(axis(box))
    ray_direction = box_invquaterion * direction(ray)
    ray_origin = box_invquaterion * (origin(ray) - origin(box))
    bounds_box = bounds(box)
    bounds(v::Float64, other::Bool=false) = @inbounds xor(signbit(v), other) ? bounds_box[2] : bounds_box[1]

    inv_ray_direction = 1.0 ./ ray_direction

    tmin = (bounds(ray_direction.x).x - ray_origin.x) * inv_ray_direction.x
    tmax = (bounds(ray_direction.x, true).x - ray_origin.x) * inv_ray_direction.x
    tymin = (bounds(ray_direction.y).y - ray_origin.y) * inv_ray_direction.y
    tymax = (bounds(ray_direction.y, true).y - ray_origin.y) * inv_ray_direction.y
    (tmin > tymax || tymin > tmax) && return nothing
    tmin = max(tmin, tymin)
    tmax = min(tmax, tymax)
    tzmin = (bounds(ray_direction.z).z - ray_origin.z) * inv_ray_direction.z
    tzmax = (bounds(ray_direction.z, true).z - ray_origin.z) * inv_ray_direction.z
    (tmin > tzmax || tzmin > tmax) && return nothing
    tmin = max(tmin, tzmin)
    tmax = min(tmax, tzmax)

    tmax > 0 || return nothing
    tmin > 0 || return distance!(minintersection, tmax)
    distance!(minintersection, tmin)
end

@approx function minintersection!(minintersection::MinIntersection, box::Box{2}, ray::Ray{2})
    hassize(box) || return nothing
    box_invquaterion = invquaterniony(axis(box))
    ray_direction = box_invquaterion * direction(ray)
    ray_origin = box_invquaterion * (origin(ray) - origin(box))
    bounds_box = bounds(box)
    bounds(v::Float64, other::Bool=false) = @inbounds xor(signbit(v), other) ? bounds_box[2] : bounds_box[1]

    inv_ray_direction = 1.0 ./ ray_direction

    tmin = (bounds(ray_direction.x).x - ray_origin.x) * inv_ray_direction.x
    tmax = (bounds(ray_direction.x, true).x - ray_origin.x) * inv_ray_direction.x
    tymin = (bounds(ray_direction.y).y - ray_origin.y) * inv_ray_direction.y
    tymax = (bounds(ray_direction.y, true).y - ray_origin.y) * inv_ray_direction.y
    (tmin > tymax || tymin > tmax) && return nothing
    tmin = max(tmin, tymin)
    tmax = min(tmax, tymax)

    0 < tmax || return nothing
    0 < tmin || return distance!(minintersection, tmax)
    distance!(minintersection, tmin)
end

function join(box1::Box{2}, box2::Box{2})
    hassize(box1) || return box2
    hassize(box2) || return box1
    b1_coordinates = coordinates(box1)
    b2_coordinates = coordinates(box2)
    box_coordinates = vcat(b1_coordinates, b2_coordinates)

    origin_x, size_x = extrema([point[1] for point in box_coordinates])
    origin_y, size_y = extrema([point[2] for point in box_coordinates])
    origin_new = SA[origin_x, origin_y]
    size_new = SA[size_x, size_y]
    Box(origin_new, size_new - origin_new, AXIS2_Y, interface(box1))
end

function join(box1::Box{3}, box2::Box{3})
    hassize(box1) || return box2
    hassize(box2) || return box1
    b1_coordinates = coordinates(box1)
    b2_coordinates = coordinates(box2)
    box_coordinates = vcat(b1_coordinates, b2_coordinates)

    origin_x, size_x = extrema([point[1] for point in box_coordinates])
    origin_y, size_y = extrema([point[2] for point in box_coordinates])
    origin_z, size_z = extrema([point[3] for point in box_coordinates])
    origin_new = SA[origin_x, origin_y, origin_z]
    size_new = SA[size_x-origin_x, size_y-origin_y, size_z-origin_z]
    Box(origin_new, size_new, AXIS3_Z, interface(box1))
end

function doesintersect(box1::Box{N}, box2::Box{N}) where {N}
    hassize(box1) || return false
    hassize(box2) || return false
    b1_edges = edges(box1)
    b2_edges = edges(box2)
    for (line_segment_1, line_segment_2) in Iterators.product(b1_edges, b2_edges)
        doesintersect(line_segment_1, line_segment_2) && return true
    end
    return false
end

function edges(box::Box{2})
    hassize(box) || return SVector{2}[]
    box_coordinates = coordinates(box)
    map(i -> (box_coordinates[i-1], box_coordinates[i]), 2:length(box_coordinates))
end

function edges(box::Box{3})
    hassize(box) || return SVector{3}[]
    box_coordinates = coordinates(box)

    edges = map(i -> (box_coordinates[i-1], box_coordinates[i]), 2:length(box_coordinates))

    missing_edges = [
        (box_coordinates[2], box_coordinates[7]),
        (box_coordinates[3], box_coordinates[8]),
        (box_coordinates[4], box_coordinates[9]),
    ]

    append!(edges, missing_edges...)
    edges
end

@approx function doesintersect(line_segment_1::NTuple{2,SVector{2,Float64}}, line_segment_2::NTuple{2,SVector{2,Float64}})
    p1, p2 = line_segment_1
    p3, p4 = line_segment_2

    A = p2 - p1
    B = p3 - p4
    C = p1 - p3

    t_numerator = B[2] * C[1] - B[1] * C[2]
    u_numerator = A[1] * C[2] - A[2]C[1]
    denominator = A[2] * B[1] - A[1] * B[2]

    if denominator > 0
        (t_numerator < 0 || u_numerator < 0) && return false
        (t_numerator > denominator || u_numerator > denominator) && return false
        return true
    end

    (t_numerator > 0 || u_numerator > 0) && return false
    (t_numerator < denominator || u_numerator < denominator) && return false

    t = t_numerator / denominator
    u = u_numerator / denominator

    p5 = p1 + t * A
    p6 = p4 + u * B

    norm(p6 - p5) == 0.0
end

@approx function doesintersect(line_segment_1::NTuple{2,SVector{3,Float64}}, line_segment_2::NTuple{2,SVector{3,Float64}})
    p11, p12 = line_segment_1
    p21, p22 = line_segment_2

    r = p12 - p11
    s = p22 - p21
    q = p11 - p21

    dotqr = dot(q, r)
    dotqs = dot(q, s)
    dotrs = dot(r, s)
    dotrr = dot(r, r)
    dotss = dot(s, s)

    numerator = dotqs * dotrs - dotqr * dotss
    denominator = dotrr * dotss - dotrs * dotrs

    if denominator < 0.0
        numerator < denominator && return false
        numerator > 0.0 && return false
        return true
    end

    numerator > denominator && return false
    numerator < 0.0 && return false

    t = numerator / denominator
    u = (dotqs + t * dotrs) / dotss

    0.0 <= u <= 1.0 || return false

    p0 = p11 + t * r
    p1 = p21 + u * s

    norm(p0 - p1) == 0.0
end

function coordinates(box::Box{2})
    hassize(box) || return SVector{2}[]
    x, y = size(box)
    points = [
        SA[0.0, 0.0],
        SA[x, 0.0],
        SA[x, y],
        SA[0.0, y],
        SA[0.0, 0.0]
    ]
    [quaterniony(axis(box)) * point + origin(box) for point in points]
end
export coordinates

function coordinates(box::Box{3})
    hassize(box) || return SVector{3}[]
    x, y, z = size(box)
    points = [
        SA[0.0, 0.0, 0.0],
        SA[x, 0.0, 0.0],
        SA[x, y, 0.0],
        SA[0.0, y, 0.0],
        SA[0.0, 0.0, 0.0],
        SA[0.0, 0.0, z],
        SA[x, 0.0, z],
        SA[x, y, z],
        SA[0.0, y, z],
        SA[0.0, 0.0, z],
    ]
    [quaternionz(axis(box)) * point + origin(box) for point in points]
end
