using Iterators
using StaticArrays
using LightXML
using WriteVTK
using Combinatorics

import Base: length

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


#
# Constructors
#
"""Curve from a series of points"""
function Curve{D, T}(points::Vector{SVector{D, T}})
    # Again, no checks!
    segments = [(i, i+1) for i in 1:(length(points)-1)] 
    Curve(points, segments)
end


"""Curve from segments defined by physical vertices"""
function Curve{D, T}(segments::Vector{Tuple{SVector{D, T}, SVector{D, T}}})
    points = Vector{SVector{D, T}}()
    index_segments = Vector{Tuple{Int, Int}}()

    count = 1
    for segment in segments
        v0, v1 = segment

        # Try finding in existing points, otherwise new
        index0 = findfirst(points, v0)
        if index0 == 0
            push!(points, v0)
            index0 = count
            count += 1
        end
        # Try finding in existing points, otherwise new
        index1 = findfirst(points, v1)
        if index1 == 0
            push!(points, v1)
            index1 = count
            count += 1
        end
        # The segment
        push!(index_segments, (index0, index1))
    end
    Curve(points, index_segments)
end

"""Load from VTU or XML file"""
function Curve(path::AbstractString)
    name, ext = splitext(path)

    ext == ".vtu" && return Curve(read_curve_vtu(path)...)

    ext == ".xml" && return Curve(read_curve_xml(path)...)

    ext == ".msh" && return Curve(read_curve_msh(path)...)   
end


"""Load from VTU"""
function read_curve_vtu(path::AbstractString)
    last(splitext(path)) == ".vtu" || throw("$(path) is not VTU file")
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
        # Finally read in
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

    # For trinagle we might be in 2d (a z=0 plane)
    nvertices == 4 && return (points, segments)
    # Not there
    norm(map(last, points), 1) > 0 && return (points, segments)
    # Chop
    points2d = Vector{SVector{2, Float64}}(npoints)
    for (row, point) in enumerate(points)
        points2d[row] = SVector(point[1], point[2])
    end

    (points2d, segments)
end

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

    read_curve_xml(path, Val{celltype})
end

    
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

    # Cell local
    edge_indices = collect(Combinatorics.combinations((1, 2, 3), 2))
    # Now we break each cell into segments    
    segments = Set{Tuple{Int, Int}}()
    for cell in imap(c -> Tuple((1+parse(Int, attribute(c, "v0")),
                                 1+parse(Int, attribute(c, "v1")),
                                 1+parse(Int, attribute(c, "v2")))),
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

    norm(map(last, points), 1) > 0 && return (points, segments)
    
    # Chop
    points2d = Vector{SVector{2, Float64}}(npoints)
    for (row, point) in enumerate(points)
        points2d[row] = SVector(point[1], point[2])
    end
    (points2d, segments)
end

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

    # Cell local
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



"""True if the curve is well defined i.e. unique points and segments"""
function is_well_defined{D, T}(c::Curve{D, T})
    # Point duplicity. Not sure if set wouldn't have been faster
    npoints = length(c.points)
    for i in 2:length(npoints)
        target = c.points[i]
        if any(target == c.points[j] for j in 1:(i-1))
            return false
        end
    end

    # Segment duplicity (orientation!) and sanity
    for i in 1:length(c.segments)
        target = c.segments[i]
        i0, i1 = target
        @assert 1 <= i0 <= npoints && 1 <= i1 <= npoints && i0 != i1
        rtarget = (i1, i0)

        if any(let
               seg_j = c.segments[j]
               target == seg_j || rtarget == seg_j
               end
               for j in 1:(i-1))
            return false
        end
    end
    # All is well
    true
end

Curve("x.xml")

#c0 = Curve("mesh2d.xml")


function write_curve{D, T}(curve::Curve{D, T}, path::AbstractString)
    name, ext = splitext(path)
    ext == ".vtu" || throw("$(path) is not an vtu file")
   
    cells = Vector{MeshCell}(length(curve.segments))
    for (i, seg) in enumerate(curve.segments)
        cells[i] = MeshCell(VTKCellTypes.VTK_LINE, [seg[1], seg[2]])
    end
    
    points = zeros((D, length(curve.points)))
    for (col, point) in enumerate(curve.points)
        for (row, value) in enumerate(point)
            points[row, col] = value
        end
    end

    vtk_save(vtk_grid(name, points, cells))
end


function koch_step(v0::SVector{2, Float64}, v1::SVector{2, Float64}, turn=:left)
    @assert v0 != v1

    v13 = v0*2/3 + v1*1/3
    v23 = v0*1/3 + v1*2/3

    t = (v1 - v0)

    shift = (turn == :left) ? SVector{2, Float64}(t[2], -t[1]) : SVector{2, Float64}(-t[2], t[1])
    
    vmid = 0.5*(v0 + v1) + shift*sqrt(3)/6
    # New points
    [v13, vmid, v23]
end


function koch(niters::Int, v0::SVector{2, Float64}, v1::SVector{2, Float64}, turn=:left)
    @assert v0 != v1
    @assert niters > 0
    # This is now pretty but perhaps faster
    points = Vector{SVector{2, Float64}}([v0, v1])
    while niters > 0
        new_points = Vector{SVector{2, Float64}}()
        for i in 1:(length(points)-1)
            v0, v1 = points[i], points[i+1]
            append!(new_points, [v0, koch_step(v0, v1, turn)...])
        end
        push!(new_points, v1)

        points = new_points
        niters -= 1
    end

    points
end


struct Segment{D, T}
    A::SVector{D, T}
    B::SVector{D, T}
end

length{D, T}(line::Segment{D, T}) = norm(line.B - line.A)

function iter_segments{D, T}(curve::Curve{D, T})
     (Segment(curve.points[i0], curve.points[i1]) for (i0, i1) in curve.segments)
end

length{D, T}(curve::Curve{D, T}) = sum(imap(length, iter_segments(curve)))

hmin{D, T}(curve::Curve{D, T}) = minimum(imap(length, iter_segments(curve)))

hmax{D, T}(curve::Curve{D, T}) = maximum(imap(length, iter_segments(curve)))

k = Curve(koch(2, SVector(0., 0), SVector(1., 1)))
println((length(k.points), length(k)))
# #write_curve(k, "test.vtu")

k = Curve(koch(4, SVector(0., 0), SVector(1., 1)))
println((length(k.points), length(k)))

Curve("mesh2d.xml")


function read_curve_msh(path::AbstractString)
    # First of all this should only work with gmsh meshes
    last(splitext(path)) == ".msh" || throw("$(path) is not an XML file")
    # For now we will only handle the following tags, others are ignored
    # FIXME: some of the intersting unhandled ones are periodic ...
    tags = ("\$MeshFormat", "\$Nodes", "\$Elements")
    # We shall always store nodes, i.e. all the coordinates used to define
    # entities of the mesh
    nodes = Vector{SVector{3, Float64}}()

    entities = Dict{Int64, Vector}()

    stream = open(path)
    # Node are always 2d but we take 2d as a z = 0 plane
    is_2d = true
    used_tag = ""
    for ln in eachline(stream)
        line = split(strip(ln), " ")
        # This is a tag. Open tags get to set used_tag, EndTag is just skipped
        if startswith(first(line), "\$")
            tag = first(line)
            if in(tag, tags)
                used_tag = tag
            end
            continue
        end
        
        # We don't do anatyhing with the line that follows line with meshformat
        if used_tag == "\$MeshFormat"
            continue
        end

        # Id Nodes tag is used
        if used_tag == "\$Nodes"
            # The array for vertices is allocated and then we move on
            if isempty(nodes)
                num_nodes = parse(Int64, first(line))
                # NOTE: GMSH stores xyz even in 2D
                resize!(nodes, num_nodes)
                continue
            end
            # The remining line fill in the coordinates
            node_index = parse(Int64, first(line))
            # Update z plane
            z = last(line)
            is_2d = (z == "0") && is_2d
            # One based indexing
            nodes[node_index] = SVector(1+parse(Float64, line[2]),
                                        1+parse(Float64, line[3]),
                                        1+parse(Float64, z))
        end

        # Lines following Element
        if used_tag == "\$Elements"
            # Are skipped with the element is point or number of element info
            if length(line) == 1 || parse(Int64, line[2]) == 15
                continue
            end

            # The rest builds gmsh elm-type
            row = map(s -> parse(Int64, s), line)
            index, elm_type, tag = row[1:3]
            elm_nodes = tuple(row[3+tag+1:end]...)

            if haskey(entities, elm_type)
                push!(entities[elm_type], elm_nodes)
            else
                entities[elm_type] = [elm_nodes]
            end
        end
    end
    close(stream)

    # At this point we have cells. Only look at those with largest dim
    celltype = maximum(keys(entities))
    @assert celltype == 2 || celltype == 4

    # Break cells into segments
    # Cell local
    range = (celltype == 2) ? (1, 2, 3) : (1, 2, 3, 4)
    edge_indices = collect(Combinatorics.combinations(range, 2))

    # Now we break each cell into segments    
    segments = Set{Tuple{Int, Int}}()
    for cell in entities[celltype]
        # Break cell to edges
        for (i0, i1) in edge_indices
            # One based indexing
            v0, v1 = cell[i0], cell[i1]
            push!(segments, (v0 < v1) ? (v0, v1) : (v1, v0))
        end
    end
    # Unique can be converted
    segments = collect(segments)

    !is_2d && return (nodes, segments)

    # Chop
    points2d = Vector{SVector{2, Float64}}(length(nodes))
    for (row, point) in enumerate(nodes)
        points2d[row] = SVector(point[1], point[2])
    end
    (points2d, segments)
end

Curve("mesh2d.msh")
