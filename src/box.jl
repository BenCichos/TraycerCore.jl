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
hassize(box::Box) = !all(iszero, size(box))

function isinside(box::Box{N}, vector::SVector{N,Float64}) where {N}
    hassize(box) || return false
    axis_vector = inv(quaterniony(axis(box))) * (vector - origin(box))

    for point in axis_vector
        @approx point < 0.0 && return false
    end

    for (point, axis_size) in zip(axis_vector, size(box))
        @approx point > axis_size && return false
    end

    return true
end

function onsurface(box::Box{N}, vector::SVector{N,Float64}) where {N}
    axis_vector = inv(quaterniony(axis(box))) * (vector - origin(box))
    for point in axis_vector
        @approx point < 0.0 && return false
    end

    for (point, axis_size) in zip(axis_vector, size(box))
        @approx point > axis_size && return false
    end

    for (point, axis_size) in zip(axis_vector, size(box))
        @approx point == 0 || point == axis_size && return true
    end

    return false
end

function normal(box::Box{2}, vector::SVector{2,Float64})
    x, y = size(box)
    axis_vector = inv(quaterniony(axis(box))) * (vector - origin(box))
    normal_x, normal_y = 0.0, 0.0
    if @approx axis_vector[1] == 0.0
        normal_x = -1.0
    elseif @approx axis_vector[1] == x
        normal_x = 1.0
    end
    if @approx axis_vector[2] == 0.0
        normal_y = -1.0
    elseif @approx axis_vector[2] == y
        normal_y = 1.0
    end

    quaterniony(axis(box)) * SA[normal_x, normal_y]
end

function normal(box::Box{3}, vector::SVector{3,Float64})
    x, y, z = size(box)
    axis_vector = inv(quaterniony(axis(box))) * (vector - origin(box))
    normal_x, normal_y, normal_z = 0.0, 0.0, 0.0
    if @approx axis_vector[1] == 0.0
        normal_x = -1.0
    elseif @approx axis_vector[1] == x
        normal_x = 1.0
    end

    if @approx axis_vector[2] == 0.0
        normal_y = -1.0
    elseif @approx axis_vector[2] == y
        normal_y = 1.0
    end

    if (@approx axis_vector[3] == 0.0)
        normal_y = -1.0
    elseif (@approx axis_vector[3] == z)
        normal_y = 1.0
    end

    quaternionz(axis(box)) * SA[normal_x, normal_y, normal_z]
end

function minintersection!(::MinIntersection, element::Box{2}, ray::Ray{2})
    hassize(element) || return nothing


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
    any(edges -> doesintersect(edges[1], edges[2]), Iterators.product(b1_edges, b2_edges)) && return true
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

function doesintersect(edge1::NTuple{2,SVector{2,Float64}}, edge2::NTuple{2,SVector{2,Float64}})
    p1, p2 = edge1
    p3, p4 = edge2

    o1 = orientation(p1, p2, p3)
    o2 = orientation(p1, p2, p4)
    o3 = orientation(p3, p4, p1)
    o4 = orientation(p3, p4, p2)

    o1 != o2 && o3 != o4 && return true

    o1 == 0 && onsegment(p1, p3, p2) && return true
    o2 == 0 && onsegment(p1, p4, p2) && return true
    o3 == 0 && onsegment(p3, p1, p4) && return true
    o4 == 0 && onsegment(p3, p2, p4) && return true

    return false
end

onsegment(p::SVector{2,Float64}, q::SVector{2,Float64}, r::SVector{2,Float64}) = @approx min(p[1], r[1]) <= q[1] <= max(p[1], r[1]) && min(p[2], r[2]) ≤ q[2] ≤ max(p[2], r[2])

function orientation(p::SVector{2,Float64}, q::SVector{2,Float64}, r::SVector{2,Float64})
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    @approx val == 0.0 && return 0
    @approx val >= 0 && return 1
    return 2
end

function doesintersect(edge1::NTuple{2,SVector{3,Float64}}, edge2::NTuple{2,SVector{3,Float64}})
    p11, p12 = edge1
    p21, p22 = edge2

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

    if @approx denominator < 0.0
        @approx numerator < denominator && return false
        @approx numerator > 0.0 && return false
        return true
    end

    @approx numerator > denominator && return false
    @approx numerator < 0.0 && return false

    t = numerator / denominator
    u = (dotqs + t * dotrs) / dotss

    @approx 0.0 <= u <= 1.0 || return false

    p0 = p11 + t * r
    p1 = p21 + u * s

    @approx norm(p0 - p1) == 0.0
end
export doesintersect

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
