using StaticArrays

const Point = SVector

"""
Curve is a collection of segments. Each segment is defined by 2 indices
corresponding to two points in the points array.
"""
struct Curve{D}
    # NOTE: it is fine to have duplicate points
    points::Vector{SVector{D}}
    # We don't shield user from putting silly indices
    segments::Vector{Tuple{Int, Int}}
end

"""Curve from a series of points"""
function Curve{D}(point::Vector{SVector{D}})
    segments = [(i, i+1) for i in 1:length(points)] 
    Curve(points, segments)
end

# 
#
#
#
#
#
#
#
x = SVector((1, 2))
y = SVector((2, 3))
points = Vector{SVector{2}}([x, y])


Curve(points)

