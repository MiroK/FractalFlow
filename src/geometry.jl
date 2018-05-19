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
include("fractal_dim.jl")

"""
Scale the curve (in place) such as the distances between points are 
changed by factor a.
"""
function scale!{D, T}(c::Curve{D, T}, a::Number)
    o = BoundingBox(c).low

    points = c.points
    for i in 1:length(c.points)
        points[i] = a*(points[i] - o) + o
    end
end


"""
The distance vector between points is scaled by diagm(a). This allows for 
sheer.
"""
function scale!{D, T}(c::Curve{D, T}, a::NTuple{D, T})
    o = BoundingBox(c).low

    X = diagm(SVector(a...))
    points = c.points
    for i in 1:length(c.points)
        points[i] = X*(points[i] - o) + o
    end
end


"""Translate! the curve by vector b"""
function translate!{D, T}(c::Curve{D, T}, b::SVector{D, T})
    points = c.points
    for i in 1:length(c.points)
        points[i] += b
    end
end


"""Scale! curve to rectangle b, i.e. the bbox of c should be b."""
function scale!{D, T}(c::Curve{D, T}, b::BoundingBox{D, T})
    target_low, target_high = b.low, b.high

    cframe = BoundingBox(c)
    c_low, c_high = cframe.low, cframe.high
    # Figure out the relative scaling
    scale = (target_high - target_low)./(c_high - c_low)

    X = diagm(scale)
    points = c.points
    for i in 1:length(c.points)
        # Scale and shift to next origin
        points[i] = X*(points[i] - c_low) + target_low
    end
end

# IDEA: fd in parallel
# TODO: tortuosity
