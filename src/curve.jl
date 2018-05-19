# Curve Constructors
"""Curve from a series of points"""
function Curve{D, T}(points::Vector{SVector{D, T}})
    # Again, no checks!
    segments = [(i, i+1) for i in 1:(length(points)-1)] 
    Curve(points, segments)
end


"""Curve from segments defined by physical vertices"""
function Curve{D, T}(segments::Vector{Tuple{SVector{D, T}, SVector{D, T}}})
    points = Dict{SVector{D, T}, Int}()
    index_segments = Set{Tuple{Int, Int}}()

    count = 1
    for segment in segments
        v0, v1 = segment

        # Try finding in existing points, otherwise new
        if v0 ∉ keys(points)
            index0 = (points[v0] = count)
            count += 1
        else
            index0 = points[v0]
        end
        
        # Try finding in existing points, otherwise new
        if v1 ∉ keys(points)
            index1 = (points[v1] = count)
            count += 1
        else
            index1 = points[v1]
        end
        edge = (index0 < index1) ? (index0, index1) : (index1, index0)
        # The segment
        push!(index_segments, edge)
    end
    points = sort(collect(keys(points)), by=k->points[k])
    index_segments = collect(index_segments)
    
    Curve(points, index_segments)
end


"""Curve from segments defined by physical vertices"""
function Curve{D, T}(segments::Vector{Segment{D, T}})
    segments = Vector{Tuple{SVector{D, T}, SVector{D, T}}}([(s.A, s.B) for s in segments])
    Curve(segments)
end


"""True if the curve is well defined i.e. unique points and segments"""
function is_consistent{D, T}(c::Curve{D, T})
    # Point duplicity. Not sure if set wouldn't have been faster
    npoints = length(c.points)
    for i in 2:npoints
        target = c.points[i]
        if any(target == c.points[j] for j in 1:(i-1))
            return false
        end
    end

    # Segment duplicity (orientation!) and sanity
    for i in 1:length(c.segments)
        target = c.segments[i]
        i0, i1 = target
        (1 <= i0 <= npoints && 1 <= i1 <= npoints && i0 != i1) || return false
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

