module FractalFlow

using StaticArrays   # SVector
using IterTools      # Imap
using LightXML       # Parse XML, VTU
using LightGraphs    #
using Combinatorics  # Combinations

import Base: length, write          # Lenght of curve, writing curve
import LightGraphs: nv, ne, Graph, degree     

export Curve, VertexFunction, EdgeFunction, is_consistent, entity_dim
export nv, ne, Graph, degree_f
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
include("topology.jl")  # Connectivity matter
# include("geometry.jl")  # Distance matters

# struct Segment{D, T}
#     A::SVector{D, T}
#     B::SVector{D, T}
# end

# length{D, T}(line::Segment{D, T}) = norm(line.B - line.A)

# function iter_segments{D, T}(curve::Curve{D, T})
#      (Segment(curve.points[i0], curve.points[i1]) for (i0, i1) in curve.segments)
# end

# length{D, T}(curve::Curve{D, T}) = sum(imap(length, iter_segments(curve)))

# hmin{D, T}(curve::Curve{D, T}) = minimum(imap(length, iter_segments(curve)))

# hmax{D, T}(curve::Curve{D, T}) = maximum(imap(length, iter_segments(curve)))

# k = Curve(koch(2, SVector(0., 0), SVector(1., 1)))
# println((length(k.points), length(k)))
# # #write_curve(k, "test.vtu")

# k = Curve(koch(4, SVector(0., 0), SVector(1., 1)))
# println((length(k.points), length(k)))

# Curve("mesh2d.xml")



# Curve("mesh2d.msh")

# # A test graph
# segments = [(SVector(0., 0.), SVector(1., 0)),
#             (SVector(1., 0.), SVector(2., 0)),
#             (SVector(0., 1.), SVector(1., 1)),
#             (SVector(1., 1.), SVector(2., 1)),
#             (SVector(0., 2.), SVector(1., 2)),
#             (SVector(1., 2.), SVector(2., 2)),
#             (SVector(0., 0.), SVector(0., 1)),
#             (SVector(0., 1.), SVector(0., 2)),
#             (SVector(1., 0.), SVector(1., 1)),
#             (SVector(1., 1.), SVector(1., 2)),
#             (SVector(2., 0.), SVector(2., 1)),
#             (SVector(2., 1.), SVector(2., 2)),
#             (SVector(0., 1.), SVector(1., 2.))]
# c1 = Curve(segments)

# using LightGraphs

# import LightGraphs: Graph

# function Graph(c::Curve)
#     g = Graph(length(c.points))
#     for seg in c.segments
#         add_edge!(g, seg...)
#     end
#     g
# end

# #####
# # Analysis
# #
# #

# function edge_vertex_connectivity(c::Curve)
#     # Make it a dict for compatibility
#     mapping = Dict{Int, Vector{Int}}()
#     for (edge, vertices) in enumerate(c.segments)
#         mapping[edge] = collect(vertices)
#     end
#     mapping
# end

# function invert{K, V}(mapping::Dict{K, Vector{V}})
#     i_mapping = Dict{V, Vector{K}}()
#     for pair in mapping
#         value = pair.first
#         for key in pair.second
#             !in(key, keys(i_mapping)) ? push!(i_mapping, Pair(key, [value])) : push!(i_mapping[key], value)
#         end
#     end
#     i_mapping
# end

# vertex_edge_connectivity(c::Curve) = invert(edge_vertex_connectivity(c))


# """A list of nodes connected to only one edge"""
# function boundary_nodes(c::Curve)
#     degrees = degree(Graph(c))
#     find(x -> x == 1, degrees)
# end

# """A list of nodes connected to more than 2 edges"""
# function branching_nodes(c::Curve)
#     degrees = degree(Graph(c))
#     find(x -> x > 2, degrees)
# end

# """Terminals [branching nodes + boundary nodes] => connected edges"""
# function terminal_map(c::Curve)
#     degrees = degree(Graph(c))
#     terminals = find(x -> x == 1 || x > 2, degrees)

#     v2e = vertex_edge_connectivity(c)
#     Dict(t => v2e[t] for t in terminals)
# end

# """Dual graph"""
# function DualGraph(g::Graph)
#     vertex_neighbors = LightGraphs.SimpleGraphs.adj(g)

#     dual_graph = Graph(ne(g))
#     for (v, neighbors) in enumerate(vertex_neighbors)
#         for n in neighbors
#             edge = (n < v) ? (n, v) : (v, n)
#             add_edge!(dual_graph, edge)
#         end
#     end
#     dual_graph
# end

# # """
# # A branch is a collection of vertices where the first and last vertices
# # are `terminal` vertices. Here we make a map of branches 
# # (first, last) => [[edges of the branch]] (where edges as given by vertex 
# # to edge map
# # """

# nnzkeys{K, V}(d::Dict{K, Vector{V}}) = filter(key -> !isempty(d[key]), keys(d))

# function branch_map(c::Curve)
#     tmap = terminal_map(c::Curve)

#     v2e = vertex_edge_connectivity(c)
#     e2v = edge_vertex_connectivity(c)

#     bmap = Dict{Tuple{Int, Int}, Vector{Vector{Int}}}()

#     while !isempty(nnzkeys(tmap))
#         # Pick som terminal
#         start = first(nnzkeys(tmap))
#         println("<< $(tmap[start])")
#         # And we will follow some edge connected to it
#         edge = pop!(tmap[start])
#         println("<<< $(tmap[start])")

#         println(start)
#         branch = [edge]

#         next_v = start
#         while true
#             println("edge $(edge) $(e2v[edge]) ")
#             v0, v1 = e2v[edge]
#             next_v = (v0 == next_v) ? v1 : v0

#             next_v ∉ nnzkeys(tmap) || break
            
#             # Next edge
#             e0, e1 = v2e[next_v]
#             println("vertex $(next_v) $(v2e[next_v])")
#             edge = (e0 in branch) ? e1 : e0
#             push!(branch, edge)
#             println("\t[$(edge)] $(next_v)")
#         end
#         # If the vertex is terminal. the path is completed by adding the edge
#         println(">>> $(branch)")
#         println("$(tmap[next_v])")
#         index = findfirst(x -> x == edge, tmap[next_v])
#         println("$(edge) $(next_v) $(tmap[next_v])")
#         index == 0 || deleteat!(tmap[next_v], index)
#         println("$(edge) $(next_v) $(tmap[next_v])")
#         key = (start < next_v) ? (start, next_v) : (next_v, start)
#         key == (start, next_v) || reverse!(branch)

#         (key in keys(bmap))? push!(bmap[key], branch) : push!(bmap, Pair(key, [branch]))
#     end
#     bmap
# end


# segments = [(SVector(0., 0.), SVector(1., 0)),
#             (SVector(1., 0.), SVector(1., -1)),
#             (SVector(1., -1.), SVector(0., -1)),
#             (SVector(0., -1.), SVector(0., 0)),
#             (SVector(0., 0.), SVector(-1., 0)),
#             (SVector(-1., 0.), SVector(-1., 1)),
#             (SVector(-1., 1.), SVector(0., 1)),
#             (SVector(0., 1.), SVector(0., 0))]

# c = Curve(segments)

# segments = [(SVector(0., 0.), SVector(1., 0)),
#             (SVector(1., 0.), SVector(1., -1)),
#             (SVector(1., -1.), SVector(0., -1)),
#             (SVector(0., -1.), SVector(0., 0)),
#             (SVector(0., 0.), SVector(-1., 0)),
#             (SVector(-1., 0.), SVector(-1., 1)),
#             (SVector(-1., 1.), SVector(0., 1)),
#             (SVector(0., 1.), SVector(0., 0)),
#             (SVector(-1., 0), SVector(0., 1)),
#             (SVector(0., -1.), SVector(1., 0))]

# cc = Curve(segments)

# # write_curve(cc, "xxx.vtu")
# using PyCall

# @pyimport networkx as nx
# @pyimport networkx.algorithms as nx_algo


# function PyGraph(g::Graph)
#     pyg = nx.Graph()
#     for e in edges(g)
#         py"$(pyg).add_edge($(e.src), $(e.dst))"
#     end
#     pyg
# end


# function number_of_loops(c::Curve)
#     bmap = branch_map(c)

#     nloops = 0
#     paths = collect(keys(bmap))
#     while any(length(v) > 1 for v in values(bmap))
#         path = findfirst(p -> length(bmap[p]) > 1, paths)
#         pop!(bmap[paths[path]])
#         nloops += 1
#     end
#     println(">>>$(nloops)")

#     # Can be represented as Graph
#     py_graph = nx.Graph()
#     for (v0, v1) in keys(bmap)
#         py"$(py_graph).add_edge($(v0), $(v1))"
#     end
    
#     nloops += size(nx_algo.cycle_basis(py_graph), 1)
#     nloops
# end


end # module
