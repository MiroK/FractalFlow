"""Euclidean norm of segment"""
length{D, T}(line::Segment{D, T}) = norm(line.B - line.A)


"""Iterator over edges of the curve as Segment entities"""
function iter_segments{D, T}(curve::Curve{D, T})
    (Segment(curve.points[i0], curve.points[i1]) for (i0, i1) in curve.segments)
end


"""Total length of the curve"""
length{D, T}(curve::Curve{D, T}) = sum(imap(length, iter_segments(curve)))


"""Smallest segment of the curve"""
min_segment{D, T}(curve::Curve{D, T}) = minimum(imap(length, iter_segments(curve)))


"""Longest segment of the curve"""
max_segment{D, T}(curve::Curve{D, T}) = maximum(imap(length, iter_segments(curve)))


"""Length of branch refered to by segments indices in curve.segments"""
function branch_length(curve::Curve, segments::Vector{Int})
    l = 0.
    for segment in segments
        # index -> indices
        v0, v1 = curve.segments[segment]
        # Vertices as physical points
        l += norm(curve.points[v1]-curve.points[v0])
    end
    l
end


"""Produce an array of branch lengths in the order they appear in bmap"""
function branch_lengths(curve::Curve,  bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    lengths = Vector{Real}(num_branches(curve))
    for key in keys(bmap)
        for branch in bmap[key]
            push!(lengths, branch_length(branch))
        end
    end
    lengths
end


"""Length of the longest branch"""
function min_branch(curve::Curve,  bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    minimum(branch_lengths(curve, bmap))
end


"""Length of the shortest branch"""
function max_branch(curve::Curve,  bmap::Dict{Tuple{Int, Int}, Vector{Vector{Int}}})
    maximum(branch_lengths(curve, bmap))
end

# Generate those which compute the map
for foo in (:branch_lengths, :min_branch, :max_branch)
    @eval $(foo)(curve::Curve) = $(foo)(curve, branch_lengths(curve))
end

######################################################################
# Box counting and fractional dimension
######################################################################

struct BoundingBox{D, T}
    low::SVector{D, T}
    high::SVector{D, T}
end

"""Bbox of segment"""
BoundingBox{D, T}(line::Segment{D, T}) = BoundingBox(SVector(Tuple(min(line.A[i], line.B[i]) for i in 1:D)),
                                                     SVector(Tuple(max(line.A[i], line.B[i]) for i in 1:D)))

                                                                    
"""Bbox of collection of segments"""
function BoundingBox{D, T}(segments::Vector{Segment{D, T}})

    low = SVector(Tuple(minimum(min(seg.A[i], seg.B[i]) for seg in segments) for i in 1:D))

    high = SVector(Tuple(maximum(max(seg.A[i], seg.B[i]) for seg in segments) for i in 1:D))

    return BoundingBox(low, high)
end


"""Bbox of a curve"""
function BoundingBox{D, T}(curve::Curve{D, T})
    segs = Vector{SVector{D, T}}(ne(curve))
    for (i, segment) in enumerate(iter_segments(curve))
        segs[i] = segment
    end
    BoundingBox(segs)
end


"""Collisions"""
# Liang-Barsky algorithm
function collides{D, T}(line::Segment{D, T}, box::BoundingBox{D, T})
    x, y = line.A
    dx = line.B[1] - x
    dy = line.B[2] - y

    p = [-dx, dx, -dy, dy]
    q = [x-box.low[1], box.high[1]-x, y-box.low[2], box.high[2]-y]

    u1 = -Inf
    u2 = Inf
    for i in 1:4
        (p[i] == 0 && q[i] < 0) && return false

        t = q[i]/p[i]

        if (p[i] < 0 && u1 < t)
            u1 = max(0, t)
        elseif (p[i] > 0 && u2 > t)
            u2 = min(1, t)
        end
    end

   (u1 > u2) && return false

    true
end

# fractional_dim_estimator
#
# fractiona_dim(c::Curve, nres)

# tortuosity
#
#

