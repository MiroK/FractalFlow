using WriteVTK   

"""Load from VTU
Obviously we want to load the line mesh. Meshes of triangles/tetrahedrons
are loaded in as a graph of edges.
"""
function read_curve_vtu(path::AbstractString)
    last(splitext(path)) == ".vtu" || throw("$(path) is not a VTU file")
    # The hierchy is doc
    #                  piece
    #                    points
    #                    cells
    xroot = root(parse_file(path))
    grid = first(child_elements(xroot))
    piece = first(child_elements(grid))

    npoints = parse(Int, attribute(piece, "NumberOfPoints"))
    
    points, cells = child_elements(piece)
    # The data is X0 Y0 Z0  X1 Y1 Z1  
    points_data = split(strip(content(points)), "  ")
    # Identify dimensionality from offets
    segments, offsets, types = child_elements(cells)

    # At this point we have one of line/tri/tet mesh. Line mesh is easiest
    # as we already have all the edges for tri and tet they need to be
    # made manually
    mesh_types = Set(split(strip(content(types)), " "))
    @assert length(mesh_types) == 1
    mesh_type = pop!(mesh_types)

    # A line mesh
    if mesh_type == "3"
        dim = parse(Int, string(first(content(offsets))))
        # Index data is i0 j0  i1 j1
        segments_data = split(strip(content(segments)), "  ")

        nsegments = parse(Int, attribute(piece, "NumberOfCells"))
p        # Finally read in
        segments = Vector{Tuple{Int, Int}}(nsegments)
        for (i, seg) in enumerate(segments_data)
            # One based indexing
            segments[i] = Tuple((1+parse(Int, v)) for v in split(seg))
        end
    
        points = Vector{SVector{dim, Float64}}(npoints)
        for (i, point) in enumerate(points_data)
            points[i] = SVector((parse(Float64, v) for v in split(point)[1:dim])...)
        end

        return (points, segments)
    end

    # A simplex mesh of triangles or tets
    @assert mesh_type == "5" || mesh_type == "10"

    nvertices = (mesh_type == "5") ? 3 : 4
    # Cell local
    edge_indices = collect(Combinatorics.combinations(collect(1:nvertices), 2))
    # Cells are made of indicices - need to split to edges
    cell_data = split(strip(content(segments)), "  ")

    # Ath this point we don't really know about the segments counts
    segments = Set{Tuple{Int, Int}}()
    for cell in imap(s -> map(x -> parse(Int, x), split(s, " ")), cell_data)
        # Break cell to edges
        for (i0, i1) in edge_indices
            # One based indexing
            v0, v1 = cell[i0]+1, cell[i1]+1
            push!(segments, (v0 < v1) ? (v0, v1) : (v1, v0))
        end
    end
    # Unique can be converted
    segments = collect(segments)

    # We load points assuming 3d
    points = Vector{SVector{3, Float64}}(npoints)
    for (i, point) in enumerate(points_data)
        points[i] = SVector((parse(Float64, v) for v in split(point))...)
    end

    # For tet the 3d assumption is correct and we're done
    nvertices == 4 && return (points, segments)
    # For triangle if it is in 3d (not z == 0 plane) 3d was right
    norm(map(last, points), 1) > 0 && return (points, segments)
    # Otherwise we project to 2d
    points2d = Vector{SVector{2, Float64}}(npoints)
    for (row, point) in enumerate(points)
        points2d[row] = SVector(point[1], point[2])
    end

    (points2d, segments)
end


"""Write curve function"""
function write(path::AbstractString, f::CurveFunction)
    name, ext = splitext(path)
    ext == ".vtu" || throw("$(path) is not a VTU file")

    curve = f.curve
    cells = Vector{MeshCell}(length(curve.segments))
    for (i, seg) in enumerate(curve.segments)
        cells[i] = MeshCell(VTKCellTypes.VTK_LINE, [seg[1], seg[2]])
    end

    points = zeros((geometrical_dim(f), length(curve.points)))
    for (col, point) in enumerate(curve.points)
        for (row, value) in enumerate(point)
            points[row, col] = value
        end
    end

    vtk_file = vtk_grid(name, points, cells)
    topological_dim(f) == 0 && vtk_point_data(vtk_file, f.data, "point_values")
    topological_dim(f) == 1 && vtk_cell_data(vtk_file, f.data, "edge_values")

    vtk_save(vtk_file)
end
