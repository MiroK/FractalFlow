using Iterators
using StaticArrays
using LightXML
using WriteVTK

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
    
    points, segments = (ext == ".vtu") ? read_curve_vtu(path) : read_curve_xml(path)
    Curve(points, segments)
end


"""Load from VTU"""
function read_curve_vtu(path::AbstractString)
    @assert splitext(path)[2] == ".vtu" || throw("$(path) is not VTU file")
    # The hierchy is doc
    #                  piece
    #                    points
    #                    cells
    xroot = root(parse_file(path))
    grid = first(child_elements(xroot))
    piece = first(child_elements(grid))

    npoints = parse(attribute(piece, "NumberOfPoints"))
    nsegments = parse(attribute(piece, "NumberOfCells"))
    
    points, cells = child_elements(piece)
    # The data is X0 Y0 Z0  X1 Y1 Z1  
    points_data = split(strip(content(points)), "  ")
    # Identify dimensionality from offets
    segments, offsets, types = child_elements(cells)
    dim = parse(string(first(content(offsets))))
    # Index data is i0 j0  i1 j1
    segments_data = split(strip(content(segments)), "  ")

    # Finally read in
    segments = Vector{Tuple{Int, Int}}(nsegments)
    for (i, seg) in enumerate(segments_data)
        # One based indexing
        segments[i] = Tuple((1+parse(v)) for v in split(seg))
    end
    
    points = Vector{SVector{dim, Float64}}(npoints)
    for (i, point) in enumerate(points_data)
        points[i] = SVector((parse(v) for v in split(point)[1:dim])...)
    end

    points, segments
end

"""Load from XML"""
function read_curve_xml(path::AbstractString)
    @assert splitext(path)[2] == ".xml" || throw("$(path) is not XML file")
    # The structure is doc
    #                    mesh
    #                     nodes
    #                     cells
    xroot = root(parse_file(path))
    mesh = first(child_elements(xroot))
    @assert attribute(mesh, "celltype") == "interval"
    
    dim = parse(attribute(mesh, "dim"))
    
    nodes, cells, _ = child_elements(mesh)

    npoints = parse(attribute(nodes, "size"))
    points = Vector{SVector{dim, Float64}}(npoints)
    # Each node is an elemn with x, y [, z] attributes
    get_coords = (dim == 2) ? (
        e -> SVector(parse(attribute(e, "x")),
                     parse(attribute(e, "y")))
    ) : (
        e -> SVector(parse(attribute(e, "x")),
                     parse(attribute(e, "y")),
                     parse(attribute(e, "z"))))
    
    for (col, elm) in enumerate(child_elements(nodes))
        points[col] = get_coords(elm)
    end

    nsegments = parse(attribute(cells, "size"))
    segments = Vector{Tuple{Int, Int}}(nsegments)
    
    for (i, cell) in enumerate(child_elements(cells))
        # One based indexing
        segments[i] = Tuple((1+parse(attribute(cell, "v0")), 1+parse(attribute(cell, "v1"))))
    end
    
    points, segments
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
        if any(target == c.points[j] || rtarget == c.points[j] for j in 1:(i-1))
            return false
        end
    end
    # All is well
    true
end

Curve("foo000000.vtu")

c0 = Curve("x.xml")


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

k = Curve(koch(1, SVector(0., 0), SVector(1., 1)))
#write_curve(k, "test.vtu")

#length(k)

           
length{D, T}(line::Tuple{SVector{D, T}, SVector{D, T}}) = norm(line[2] - line[1])

function iter_segments{D, T}(curve::Curve{D, T})
    ((curve.points[i0], curve.points[i1]) for (i0, i1) in curve.segments)
end

length{D, T}(curve::Curve{D, T}) = sum(imap(length, iter_segments(curve)))

hmin{D, T}(curve::Curve{D, T}) = minimum(imap(length, iter_segments(curve)))

hmax{D, T}(curve::Curve{D, T}) = maximum(imap(length, iter_segments(curve)))

                       
