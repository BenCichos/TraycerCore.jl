mutable struct MinIntersection
    distance::Float64
    index::Int
    update_index::Bool
    intersection::Intersection

    MinIntersection() = new(Inf, -1, false, Intersection())
end

distance(minintersection::MinIntersection) = minintersection.distance
index(minintersection::MinIntersection) = minintersection.index
intersection(minintersection::MinIntersection) = minintersection.intersection

function distance!(minintersection::MinIntersection, newdistance::Float64)
    distance(minintersection) <= newdistance && return nothing
    minintersection.distance = newdistance
    minintersection.update_index = true
    return nothing
end

function update_index!(minintersection::MinIntersection, newindex::Int)
    minintersection.update_index || return nothing
    minintersection.index = newindex
    minintersection.update_index = false
    return nothing
end

reset!(minintersection::MinIntersection) = (minintersection.distance = Inf; minintersection.index = -1; minintersection.update_index = false; return nothing)

minintersection!(::MinIntersection, element::AbstractOpticalElement{N}, ray::Ray{N}) where {N} = throw(ErrorException("No method matching minintersection!(::MinIntersection, ::$(typeof(element)), ::$(typeof(ray)))"))
minintersection!(minint::MinIntersection, c::CEP, r::Ray{N}) where {N,CEP<:CompoundElementPrimitive{N}} = minintersection!(minint, primitive_element(c), r)

export MinIntersection
export distance!, update_index!, reset!, minintersection!
