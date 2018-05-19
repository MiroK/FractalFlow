# We want to generate a geo file for curve surrounded by a rectangle
# defined by shift to the bounding box. If the box is tight some of the
# the curve's points will be on the boundary and we need to help gmsh
# to correctly draw the rectanngle

"""Distance between interior point p and box""" 
function distance{D, T}(b::BoundingBox{D, T}, p::SVector{D, T})
    lx, hx = b.low[1], b.high[1]
    (lx <= p[1] <= hx) || throw("Point not inside")

    ly, hy = b.low[2], b.high[2]
    (ly <= p[2] <= hy) || throw("Point not inside")

    min(p[1] - lx, hx - p[1], p[2] - ly, hy - p[2])
end

"""Order points such that when visited one after other a loop is formed"""
function loop{T}(points::Vector{SVector{2, T}})
    # The idea here is to order according to angle - think polar coordinates
    x = map(first, points)
    x0 = mean(x)

    y = map(last, points)
    y0 = mean(y)

    angles = atan2.(x, y)
    loop = sortperm(angles)
    # Adding the last guy to convenienly makes lines with zip
    push!(loop, first(loop))
end


"""Generate a geo file for code contained in a padded bounding box"""
function gmsh_code{T}(c::Curve{2, T},
                      xpadding::Tuple{Number, Number},
                      ypadding::Tuple{Number, Number},
                      path::AbstractString)
    # Nongetative
    @assert first(xpadding) >= zero(T) && last(xpadding) >= zero(T)
    @assert first(ypadding) >= zero(T) && last(ypadding) >= zero(T)
    
    box = BoundingBox(c)
    low, high = box.low, box.high
    # Enclosing frame for the curve; shift to bbox
    frame_pts = [low - SVector(first(xpadding), first(ypadding)),
                 SVector(high[1], low[2]) + SVector(last(xpadding), -first(ypadding)),
                 box.high + SVector(last(xpadding), last(ypadding)),
                 SVector(low[1], high[2]) + SVector(-first(xpadding), last(ypadding))]

    # With points on a boundary we must think how to define the boundaing
    # rectangle in gmsh
    if !(all(xpadding .> 0) && all(ypadding .> 0))
        box = BoundingBox(frame_pts[1], frame_pts[3])
        # Now let's get all the points that are on the boundary
        d = p -> distance(box, p)
        distances = map(d, c.points)
        indices = sortperm(distances)
        # Visiting indices (ordered by distance)
        bdry_indices = Vector{Int}()
        for i in indices
            # Since we have ordering the first guys is not on bondary
            # is also the last one checked
            distances[i] > 0 && break
            # On boundary
            push!(bdry_indices, i)
        end


        bdry_pts = c.points[bdry_indices]
        # It might be that the frame points are already among the boundary ones
        # Then we wat to remove them
        for i in 1:length(frame_pts)
            last(frame_pts) ∈ bdry_pts && pop!(frame_pts)
        end

        # Using bdry and frame points we wamt to draw a loop
        loop_pts = Vector{SVector{2, T}}([bdry_pts..., frame_pts...])
        # Get the loop refered to by coordinates of loop_pts (local)
        loop_idx = loop(loop_pts)
        nbdry, nall = length(bdry_pts), length(c.points)
        # Translate to coords of curve points
        # If local refers to bdry point do a lookup for the global
        # otherwise the numbering is started after the curve points
        loop_idx = [(p <= nbdry) ? bdry_indices[p] : (nall + p - nbdry)
                    for p in loop_idx]
        # And the loop
        bdry_lines = collect(zip(loop_idx[1:end-1], loop_idx[2:end]))
    else
        n = length(c.points)
        bdry_lines = [(1+n, 2+n), (2+n, 3+n), (3+n, 4+n), (4+n, 1+n)]
    end
    # Done with prep work"""
    gmsh_code(c, frame_pts, bdry_lines, path)
end

"""Write the gmsh code"""
function gmsh_code{T}(c::Curve,
                      bdry_pts::Vector{SVector{2, T}},
                      bdry_lines::Vector{Tuple{Int, Int}},
                      path::AbstractString)

    f = open(path, "w")
    # Allow for different sizes on bdry and curve
    # NOTE: call gmsh -setnumber size VALUE0 -setnumber fsize VALUE1 ...
    write(f, "DefineConstant[size = {1}];\n")
    write(f, "DefineConstant[fsize = {1}];\n")

    # Write curve first
    for (index, p) in enumerate(c.points)
        write(f, "Point($(index)) = {$(p[1]), $(p[2]), 0, size};\n")
    end

    for (index, line) in enumerate(c.segments)
        write(f, "Line($(index)) = {$(line[1]), $(line[2])};\n")
    end
    # Mark crack as physical line 1
    crack = join(map(string, 1:length(c.segments)), ", ")
    write(f, "crack[] = {$(crack)};\n")
    write(f, "Physical Line(1) = {crack[]};\n")

    # Write bdry, points/lines start their enumeration after = consistent
    # with how bdry points are numbered in bdry lines
    n = length(c.points)
    for (i, p) in enumerate(bdry_pts)
        index = i + n
        write(f, "Point($(index)) = {$(p[1]), $(p[2]), 0, fsize};\n")
    end

    n = length(c.segments)
    for (i, line) in enumerate(bdry_lines)
        index = i + n
        write(f, "Line($(index)) = {$(line[1]), $(line[2])};\n")
    end

    loop = join(map(string, (n+1):(n+length(bdry_lines))), ", ")
    write(f, "Curve Loop(1) = {$(loop)};\n")
    write(f, "Plane Surface(1) = {1};\n")
    write(f, "Physical Surface(1) = {1};\n")
    # The magic
    write(f, "Line{crack[]} In Surface{1};\n")

    close(f)
end
