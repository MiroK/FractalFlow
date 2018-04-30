module FractalFlow

using StaticArrays   # SVector
using IterTools      # Imap
using LightXML       # Parse XML, VTU
using LightGraphs    #
using Combinatorics  # Combinations
using DataStructures # Deque
using Polynomials


import Base: length, write, size, split, start, next, done
import LightGraphs: nv, ne, Graph, degree     
# Topology
export Curve, VertexFunction, EdgeFunction, is_consistent, entity_dim
export nv, ne, Graph, degree_f, num_loops, branch_map, num_branches, num_colors, color_f
# Geometry
export FractalDim, fractal_dim, estimate_fd, min_branch, max_branch
export min_segment, max_segment
# Generation, fractals
export koch_t_curve, koch_q1_curve, koch_q2_curve
# Generation, random
# IO
export write

# -------------------------------------------------------------------

"""
Curve is a collection of segments. Each segment is defined by 2 indices
corresponding to two points in the points array.
"""
struct Curve{D, T}
    # NOTE: The inner preconditioner assumes sensible data that is points
    # in the vector are unique and the segments are unique as well. We do
    # no check if this is true here
    points::Vector{SVector{D, T}}
    segments::Vector{Tuple{Int, Int}}
end
geometrical_dim{D, T}(::Curve{D, T}) = D


# Constructors
include("curve.jl")

abstract type CurveFunction end

geometrical_dim(f::CurveFunction) = geometrical_dim(f.curve)


"""Values in nodes"""
struct VertexFunction{R <: Real} <: CurveFunction
    curve::Curve
    data::Vector{R}

    function VertexFunction{R}(c::Curve, d::Vector{R}) where R <: Real
        @assert length(c.points) == length(d)
        new(c, d)
    end
end
VertexFunction(c::Curve, d::Vector) = VertexFunction{eltype(d)}(c, d)

topological_dim(::VertexFunction) = 0


"""Values in segments"""
struct EdgeFunction{R <: Real} <: CurveFunction
    curve::Curve
    data::Vector{R}

    function EdgeFunction{R}(c::Curve, d::Vector{R}) where R <: Real
        @assert length(c.segments) == length(d)
        new(c, d)
    end
end
EdgeFunction(c::Curve, d::Vector) = EdgeFunction{eltype(d)}(c, d) 

topological_dim(::EdgeFunction) = 1

# IO - let curves be 'readable' from files
include("vtk_io.jl")
include("xml_io.jl")
include("msh_io.jl")

"""Load from VTU or XML file"""
function Curve(path::AbstractString)
    name, ext = splitext(path)

    ext == ".vtu" && return Curve(read_curve_vtu(path)...)
    
    ext == ".xml" && return Curve(read_curve_xml(path)...)

    ext == ".msh" && return Curve(read_curve_msh(path)...)   
end

# Various properties of curves
include("topology.jl")  # Connectivity

"""|AB|"""
struct Segment{D, T}
    A::SVector{D, T}
    B::SVector{D, T}
end

include("geometry.jl")  # Distance, fractal properties


# Generation
include("generation.jl")  # Fractals and random curves

end # module
