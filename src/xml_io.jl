"""Load from XML"""
function read_curve_xml(path::AbstractString)
    last(splitext(path)) == ".xml" || throw("$(path) is not XML file")
    # The structure is doc
    #                    mesh
    #                     nodes
    #                     cells
    xroot = root(parse_file(path))
    mesh = first(child_elements(xroot))

    celltype = Symbol(attribute(mesh, "celltype"))
    # We dispatch on the cell type: interval or else
    read_curve_xml(path, Val{celltype})
end

"""Load line mesh as points and segments(by indices)"""
function read_curve_xml(path::AbstractString, cell_type::Type{Val{:interval}})    
    last(splitext(path)) == ".xml" || throw("$(path) is not XML file")
    # The structure is doc
    #                    mesh
    #                     nodes
    #                     cells
    xroot = root(parse_file(path))
    mesh = first(child_elements(xroot))

    dim = parse(Int, attribute(mesh, "dim"))
    
    nodes, cells = child_elements(mesh)

    npoints = parse(Int, attribute(nodes, "size"))
    points = Vector{SVector{dim, Float64}}(npoints)
    # Each node is an elemn with x, y [, z] attributes
    get_coords = (dim == 2) ? (
        e -> SVector(parse(Float64, attribute(e, "x")),
                     parse(Float64, attribute(e, "y")))
    ) : (
        e -> SVector(parse(Float64, attribute(e, "x")),
                     parse(Float64, attribute(e, "y")),
                     parse(Float64, attribute(e, "z"))))
    
    for (col, elm) in enumerate(child_elements(nodes))
        points[col] = get_coords(elm)
    end

    nsegments = parse(Int, attribute(cells, "size"))
    segments = Vector{Tuple{Int, Int}}(nsegments)
    
    for (i, cell) in enumerate(child_elements(cells))
        # One based indexing
        segments[i] = Tuple((1+parse(Int, attribute(cell, "v0")), 1+parse(Int, attribute(cell, "v1"))))
    end
    
    (points, segments)
end

"""Turn a triangle mesh into the edge graph mesh"""
function read_curve_xml(path::AbstractString, cell_type::Type{Val{:triangle}})    
    last(splitext(path)) == ".xml" || throw("$(path) is not XML file")
    # The structure is doc
    #                    mesh
    #                     nodes
    #                     cells
    xroot = root(parse_file(path))
    mesh = first(child_elements(xroot))

    nodes, cells = child_elements(mesh)

    npoints = parse(Int, attribute(nodes, "size"))
    points = Vector{SVector{3, Float64}}(npoints)    
    # Each node is an elemn with x, y z attributes    
    for (col, elm) in enumerate(child_elements(nodes))
        points[col] = SVector(parse(Float64, attribute(elm, "x")),
                              parse(Float64, attribute(elm, "y")),
                              parse(Float64, attribute(elm, "z")))
    end

    # Local indides of triangle
    edge_indices = collect(Combinatorics.combinations((1, 2, 3), 2))
    # Now we break each cell into segments    
    segments = Set{Tuple{Int, Int}}()
    for cell in imap(c -> Tuple((1+parse(Int, attribute(c, "v0")),
                                 1+parse(Int, attribute(c, "v1")),
                                 1+parse(Int, attribute(c, "v2")))),
                     child_elements(cells))
        # All the edges
        for (i0, i1) in edge_indices
            # One based indexing
            v0, v1 = cell[i0], cell[i1]   # Global
            push!(segments, (v0 < v1) ? (v0, v1) : (v1, v0))
        end
    end
    # Unique can be converted
    segments = collect(segments)
    # If we are not in the z plane
    norm(map(last, points), 1) > 0 && return (points, segments)
    
    # Project to 2d
    points2d = Vector{SVector{2, Float64}}(npoints)
    for (row, point) in enumerate(points)
        points2d[row] = SVector(point[1], point[2])
    end
    (points2d, segments)
end

"""Turn a tetrahedron mesh into the edge graph mesh"""
function read_curve_xml(path::AbstractString, cell_type::Type{Val{:tetrahedron}})    
    last(splitext(path)) == ".xml" || throw("$(path) is not XML file")
    # The structure is doc
    #                    mesh
    #                     nodes
    #                     cells
    xroot = root(parse_file(path))
    mesh = first(child_elements(xroot))

    nodes, cells = child_elements(mesh)

    npoints = parse(Int, attribute(nodes, "size"))
    points = Vector{SVector{3, Float64}}(npoints)    
    # Each node is an elemn with x, y z attributes    
    for (col, elm) in enumerate(child_elements(nodes))
        points[col] = SVector(parse(Float64, attribute(elm, "x")),
                              parse(Float64, attribute(elm, "y")),
                              parse(Float64, attribute(elm, "z")))
    end

    # Cell local vertex indices
    edge_indices = collect(Combinatorics.combinations((1, 2, 3, 4), 2))
    # Now we break each cell into segments    
    segments = Set{Tuple{Int, Int}}()
    for cell in imap(c -> Tuple((1+parse(Int, attribute(c, "v0")),
                                 1+parse(Int, attribute(c, "v1")),
                                 1+parse(Int, attribute(c, "v2")),
                                 1+parse(Int, attribute(c, "v3")))),
                     child_elements(cells))
        # Break cell to edges
        for (i0, i1) in edge_indices
            # One based indexing
            v0, v1 = cell[i0], cell[i1]
            push!(segments, (v0 < v1) ? (v0, v1) : (v1, v0))
        end
    end
    # Unique can be converted
    segments = collect(segments)

    (points, segments)
end
