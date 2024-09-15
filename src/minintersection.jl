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
    if distance(minintersection) > newdistance
        minintersection.distance = newdistance
        minintersection.update_index = true
    end
    return nothing
end

function update_index!(minintersection::MinIntersection, newindex::Int)
    if minintersection.update_index
        minintersection.index = newindex
        minintersection.update_index = false
    end
    return nothing
end

reset!(minintersection::MinIntersection) = (minintersection.distance = Inf; minintersection.index = -1; minintersection.update_index = false; return nothing)

function minintersection!(::MinIntersection, element::AbstractOpticalElement{N}, ray::Ray{N}) where {N}
    throw(ErrorException("NotImplementedError: no method matching minintersection!(::MinIntersection, ::$(typeof(element)), ::$(typeof(ray))) has been implemented"))
    return nothing
end

minintersection!(minintersection::MinIntersection, compound_element_primitive::CEP, ray::Ray{N}) where {N,CEP<:CompoundElementPrimitive{N}} = minintersection!(minintersection, compound_element_primitive.primitive_element, ray)
