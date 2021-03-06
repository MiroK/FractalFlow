"""
Get the curve data from MSH file by breaking triangles and test into
segments.
"""
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
    # Assert simplices tris or tets
    @assert celltype == 2 || celltype == 4

    # Break cells into segments
    # Cell local index numbering
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
    # For tets or embedded tris we're done
    !is_2d && return (nodes, segments)

    # Project to 2d
    points2d = Vector{SVector{2, Float64}}(length(nodes))
    for (row, point) in enumerate(nodes)
        points2d[row] = SVector(point[1], point[2])
    end
    (points2d, segments)
end
