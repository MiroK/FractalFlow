module FractalFlow

using StaticArrays

export Point, Curve

const Point = SVector

# Collection of segments
include("curve.jl")

end # module
