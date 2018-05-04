using PyCall

@pyimport networkx as nx
@pyimport networkx.algorithms as nx_algo


"""Number of vertices in the curve"""
nv(c::Curve) = length(c.points)

"""Number of edges in the curve"""
ne(c::Curve) = length(c.segments)

"""Map edge indices to the 2 vertex indices connected to it"""
function edge_vertex_connectivity(c::Curve)
    # Make it a dict for compatibility
    mapping = Dict{Int, Vector{Int}}()
    for (edge, vertices) in enumerate(c.segments)
        mapping[edge] = collect(vertices)
    end
    mapping
end

# FIXME
"""
Compute the inverse mapping of when single K maps to several Vs. So 
one V maps to several Ks.
"""
function invert{K, V}(mapping::Dict{K, Vector{V}})
    i_mapping = Dict{V, Set{K}}()
    for pair in mapping
        value = pair.first
        for key in pair.second
            haskey(i_mapping, key) ? push!(i_mapping[key], value) : push!(i_mapping, Pair(key, Set(value)))
        end
    end
    Dict(k => collect(v) for (k, v) in i_mapping)
end


"""Map of vertex indices to at least one edge that is connected to it"""
vertex_edge_connectivity(c::Curve) = invert(edge_vertex_connectivity(c))

# For several properties it will be useful to rely on representation of
# the curve as graph
function Graph(c::Curve)
    g = Graph(nv(c))
    for edge in c.segments
        add_edge!(g, edge)
    end
    # Graph is not MultiGraph!
    @assert ne(c) == ne(g)
    g
end


"""A list of nodes [indices in vertices(graph)] connected to only one edge"""
function boundary_vertices(g::Graph)
    degrees = degree(g)
    find(x -> x == 1, degrees)
end

boundary_vertices(c::Curve) = boundary_vertices(Graph(c))


"""A list of nodes [indices in vertice(graph)] connected to more than 2 edges"""
function branching_vertices(g::Graph)
    degrees = degree(g)
    find(x -> x > 2, degrees)
end

branching_vertices(c::Curve) = branching_vertices(Graph(c))

"""An array holding for each vertex[index] a number of edges connected to it"""
degree(c::Curve) = degree(Graph(c))

"""Represent degree of the curve as VertexFunction"""
degree_f(c::Curve) = VertexFunction(c, degree(c))


"""
We are interested in branching of the curve. A branch is a collection 
of edges bounded by 2 vertices which are terminal, meaning that the vertex 
is either a branching vertex or a boundary vertex. Here we pair the terminal 
vertices with the edges connected to them.
"""
function terminal_map(c::Curve)
    degrees = degree(c)
    # Avoid building the graph twice
    terminals = find(x -> x == 1 || x > 2, degrees)

    v2e = vertex_edge_connectivity(c)
    Dict(t => v2e[t] for t in terminals)
end


"""
A branch is a collection of vertices where the first and last vertices
are `terminal` vertices. Here we make a map of branches 
(first, last) => [[edges of the branch]] (where edges as given by vertex 
to edge map. NOTE the value datatype is Vector{Vector{Int}} because in 
general there can be multiple branches between the vertices. In such case 
the reduced curve (each branch is an edge) cannot be represented as a Graph
(but MultiGraph is okay)
"""
function branch_map(c::Curve)
    tmap = terminal_map(c::Curve)

    v2e = vertex_edge_connectivity(c)
    e2v = edge_vertex_connectivity(c)

    terminals = Dict(t => length(tmap[t]) for t in keys(tmap))
    # The idea is to pich a terminal, follow some edge connected to 
    # it until we hit a terminal. We keep removing the starting/ending
    # edges so that the at some point nnzkeys is empty
    bmap = Dict{Tuple{Int, Int}, Vector{Vector{Int}}}()
    while !isempty(terminals)
        # Pick som terminal
        start = first(keys(terminals))
        # And we will follow some edge connected to it; leaving no trail
        edge = pop!(tmap[start])
        terminals[start] -= 1
        terminals[start] == 0 && delete!(terminals, start)  
        
        branch = [edge]
        next_v = start
        while true
            v0, v1 = e2v[edge]
            # Pick a different vertex then which lead us here
            next_v = (v0 == next_v) ? v1 : v0
            # Terminal vertex?
            next_v ∈ keys(terminals) && break
            
            # Otherwise this is a vertex with 2 edges. Next edge is different
            # Then current
            e0, e1 = v2e[next_v]
            edge = (e0 == edge) ? e1 : e0
            push!(branch, edge)
        end
        # If the vertex is terminal remove the edge from the end so that
        # we don't follow the path back
        deleteat!(tmap[next_v], findfirst(x -> x == edge, tmap[next_v]))
        # This reduces connections of the terminal
        terminals[next_v] -= 1
        terminals[next_v] == 0 && delete!(terminals, next_v)

        # We make a rule that key is an ordered pair 
        key = (start < next_v) ? (start, next_v) : (next_v, start)
        # And adjust the path if we flipped
        key != (start, next_v) && reverse!(branch)

        (key in keys(bmap))? push!(bmap[key], branch) : push!(bmap, Pair(key, [branch]))
    end
    bmap
end


"""Branch count"""
function num_branches(c::Curve, bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    sum((length(v) for v in values(bmap)))
end

num_branches(c::Curve) = num_branches(c, branch_map(c))


"""
Loops (as collection of edges) in the curve (with precomputed branch map)
NOTE: We produce independent basis not the `minimal` basis by which I mean
that

 ___
| /|  
|/_| is not decomposed into 2 triangles but maybe one tri + a square.
"""
function loop_basis!(c::Curve, bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    loops = Vector{Vector{Int}}()

    # Self loops (start, start) can be added right-away
    paths = collect(keys(bmap))
    closed = find(p -> p[1] == p[2], paths)
    for c in reverse!(closed)
        for loop in bmap[paths[c]]
            push!(loops, loop)
        end
        delete!(bmap, paths[c])
    end

    # To reduce to graph we build loops from multi edges
    paths = collect(keys(bmap))
    while any(length(v) > 1 for v in values(bmap))
        # Get a multiedge
        path = paths[findfirst(p -> length(bmap[p]) > 1, paths)]
        # This could be a pop, but I will remove the one which is longest
        sort!(bmap[path], by=length)
        there = pop!(bmap[path])
        # Each removed edge formed a loop with some of the remaining
        # This is an okay question
        back = bmap[path][end]
        # A loop there - back
        push!(loops, append!(there, back))
    end
    
    # NOTE: LightGraphs does not have cycle basis so we reach out to python
    py_graph = nx.Graph()
    for (v0, v1) in keys(bmap)
        py"$(py_graph).add_edge($(v0), $(v1))"
    end

    # Here in each row is 
    cycle_basis = nx_algo.cycle_basis(py_graph)
    for row in 1:size(cycle_basis, 1)
        vertices = cycle_basis[row, :]

        # [1, 2, 3, 4] = (1, 2), (2, 3), (3, 4), (4, 1)
        loop = Vector{Int}()
        for (start, stop) in zip(vertices, vcat(vertices[2:end], vertices[1]))
            # We have sorted keys
            path = (start < stop) ? (start, stop) : (stop, start)
            @assert length(bmap[path]) == 1
            append!(loop, first(bmap[path]))
        end
        push!(loops, loop)
    end
    loops
end

"""Loops in the curve (without precomputed branch map)"""
loop_basis(c::Curve) = loop_basis!(c, branch_map(c))


"""Number of loops"""
function num_loops(c::Curve, bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    # The number of loops in the (vertex)graph is the same as the number
    # of zeros in the spectrum of the laplacian matrix of the dual(edge) graph.
    # This definition is potentialy very expensive to use so the here is
    # to first buld a reduced graph where each branch in the orginal
    # graph is a signel edges => smaller graph. However, in the reduction
    # process the graph might become a MultiGraph so we first reduce
    # to Graph case

    #  __
    # |  |
    # |__|
    nloops = 0
    paths = collect(keys(bmap))
    # The case where the situation is obcious from keys
    closed = find(p -> p[1] == p[2], paths)
    for c in reverse!(closed)
        nloops += length(bmap[paths[c]])
        deleteat!(paths, c)
    end
    
    #  __ B
    # | /|
    # |/_|
    # A
    multi_edges = filter(p -> length(bmap[p]) > 1, paths)
    # Count multiplicity of multi edges
    nloops += isempty(multi_edges) ? 0 : sum(length(bmap[edge])-1 for edge in multi_edges) 

    # NOTE: LightGraphs does not have cycle basis so we reach out to python
    py_graph = nx.Graph()
    for (v0, v1) in paths
        py"$(py_graph).add_edge($(v0), $(v1))"
    end
    # Here in each row is one cycle basis
    nloops += size(nx_algo.cycle_basis(py_graph), 1)
    nloops
end

num_loops(c::Curve) = num_loops(c, branch_map(c))

"""
Graph from the mesh has as edges the line elements. In dual graph the 
edges are vertices. Zeros of the Laplace(DualGraph(curve)) count the 
loops
"""
function DualGraph(g::Graph)
    vertex_neighbors = LightGraphs.SimpleGraphs.adj(g)

    dual_graph = Graph(ne(g))
    for (v, neighbors) in enumerate(vertex_neighbors)
        for n in neighbors
            add_edge!(dual_graph, (n, v))
        end
    end
    dual_graph
end

DualGraph(c::Curve) = DualGraph(Graph(c))

"""
If the graph is split into branches how many colors are needed to 
colorthem so that no 2 neighboring branches have the same color.
"""
function num_colors(c::Curve, bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    # Build a reduced representation where each branch is a vertex
    # a graph. We code the branch as bmap[first(code)][last(code)]
    codes = Vector{Tuple{Tuple{Int, Int}, Int}}()
    for key in keys(bmap)
        for vertex in 1:length(bmap[key])
            push!(codes, (key, vertex))
        end
    end

    # Edges
    g = Graph(length(codes))
    for (v0, code0) in enumerate(codes)
        e0, e1 = first(code0)  # Vertices of the branch
        for (v1, code1) in enumerate(codes[(v0+1):end])
            if (e0 in first(code1) || e1 in first(code1))
                add_edge!(g, (v0, v0+v1))
            end 
        end
    end
    
    colors = greedy_color(g)
    colors.num_colors
end

num_colors(c::Curve) = num_colors(c, branch_map(c))

"""
If the graph is split into branches how many colors are needed to 
colorthem so that no 2 neighboring branches have the same color. 
Produce a color function for the graph (EdgeFunction).
"""
function color_f(curve::Curve, bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    # Build a reduced representation where each branch is a vertex
    # a graph. We code the branch as bmap[first(code)][last(code)]
    codes = Vector{Tuple{Tuple{Int, Int}, Int}}()
    for key in keys(bmap)
        for vertex in 1:length(bmap[key])
            push!(codes, (key, vertex))
        end
    end

    # Edges
    g = Graph(length(codes))
    for (v0, code0) in enumerate(codes)
        e0, e1 = first(code0)  # Vertices of the branch
        for (v1, code1) in enumerate(codes[(v0+1):end])
            if (e0 in first(code1) || e1 in first(code1))
                add_edge!(g, (v0, v0+v1))
            end 
        end
    end
    color_map = greedy_color(g)

    color_data = Vector{Int}(ne(curve))
    for (branch, color) in enumerate(color_map.colors)
        key, index = codes[branch]
        for edge in bmap[key][index]
            color_data[edge] = color
        end
    end
    EdgeFunction(curve, color_data)
end

color_f(c::Curve) = color_f(c, branch_map(c))


#using StaticArrays 
segments = [(SVector(0., 0.), SVector(1., 0)),
            (SVector(1., 0.), SVector(1., -1)),
            (SVector(1., -1.), SVector(0., -1)),
            (SVector(0., -1.), SVector(0., 0)),
            (SVector(0., 0.), SVector(-1., 0)),
            (SVector(-1., 0.), SVector(-1., 1)),
            (SVector(-1., 1.), SVector(0., 1)),
            (SVector(0., 1.), SVector(0., 0))]
c = Curve(segments)


segments = [(SVector(0., 0.), SVector(1., 0)),
            (SVector(1., 0.), SVector(1., -1)),
            (SVector(1., -1.), SVector(0., -1)),
            (SVector(0., -1.), SVector(0., 0)),
            (SVector(0., 0.), SVector(-1., 0)),
            (SVector(-1., 0.), SVector(-1., 1)),
            (SVector(-1., 1.), SVector(0., 1)),
            (SVector(0., 1.), SVector(0., 0)),
            (SVector(-1., 0), SVector(0., 1)),
            (SVector(0., -1.), SVector(1., 0))]
cc = Curve(segments)

segments = [(SVector(0., 0.), SVector(1., 0)),
            (SVector(1., 0.), SVector(1., -1)),
            (SVector(1., -1.), SVector(0., -1)),
            (SVector(0., -1.), SVector(0., 0)),
            (SVector(0., 0.), SVector(-1., 0)),
            (SVector(-1., 0.), SVector(-1., 1)),
            (SVector(-1., 1.), SVector(0., 1)),
            (SVector(0., 1.), SVector(0., 0)),
            (SVector(-1., 0), SVector(0., 1)),
            (SVector(0., -1.), SVector(1., 0)),
            (SVector(1., -1.), SVector(1., -2)),
            (SVector(1., -2.), SVector(2., -2)),
            (SVector(2., -2.), SVector(2., -1)),
            (SVector(2., -1.), SVector(1., -1))]
ccc = Curve(segments)

segments = [(SVector(0., 0.), SVector(1., 0)),
            (SVector(1., 0.), SVector(1., 1)),
            (SVector(1., 1.), SVector(0., 1)),
            (SVector(0., 1.), SVector(0., 0)),
            (SVector(0., 0.), SVector(1., 1.))]
cccc = Curve(segments)
