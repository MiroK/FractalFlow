# The goal is to estimate the fractal dimension D of a curve by box counting
# algorith as N = eps^{-D}, where N is the number of boxes intersected
# by the curve at this level of refinement eps.

"""A hyper prism"""
struct BoundingBox{D, T}
    low::SVector{D, T}
    high::SVector{D, T}
end


"""Physical size of the bounding box; D vector of sizes"""
size{D, T}(box::BoundingBox{D, T}) = box.high - box.low


"""Length of bounding box in the i-th axis"""
size{D, T}(box::BoundingBox{D, T}, i::Int) = size(box)[i]


"""One refinement of a 2d box results in 4 boxes"""
function split{T}(box::BoundingBox{2, T})
    X, Y = box.low, box.high
    M = 0.5(X + Y)
    [BoundingBox(X, M),
     BoundingBox(SVector{2, T}(X[1], M[2]), SVector{2, T}(M[1], Y[2])),
     BoundingBox(SVector{2, T}(M[1], X[2]), SVector{2, T}(Y[1], M[2])),
     BoundingBox(M, Y)]
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
    segs = Vector{Segment{D, T}}(ne(curve))
    for (i, segment) in enumerate(iter_segments(curve))
        segs[i] = segment
    end
    BoundingBox(segs)
end


"""Does line collide with a box?"""
# NOTE: Liang-Barsky algorithm
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

"""
Box counter holds after each refinement a pair of curve and box which 
might intersect it. The curve is a collection if indices which refer 
to segments of the unsplit curve.
"""
struct BoxCounter{D, T}
    segments::Vector{Segment{D, T}}
    seg_box::Deque{Pair{Vector{Int}, BoundingBox{D, T}}}
    stop_length::Float64
end


"""Box size at the current level"""
size{D, T}(bc::BoxCounter{D, T}) = size(first(bc.seg_box).second)


"""Construct BoxCounter for the curve"""
function BoxCounter{D, T}(curve::Curve{D, T})
    segs = Vector{Segment{D, T}}(ne(curve))
    for (i, segment) in enumerate(iter_segments(curve))
        segs[i] = segment
    end
    box = BoundingBox(segs)
    # Make the box slightly larger so that we don't hit corners
    X, Y = box.low, box.high
    center = (X+Y)/2.
    
    # Increaze the box a bit in an attempt/delay point isects
    # NOTE: as a result of this each run gives a slighly different
    # number - ideally avarage several runs
    shift1 = SVector(((Y - center).*(1. + rand(D)))...)
    shift2 = SVector(((Y - center).*(1. + rand(D)))...)
    
    box = BoundingBox(center - shift1, center + shift2)

    seg_box = Deque{Pair{Vector{Int}, BoundingBox{D, T}}}()
    push!(seg_box, Pair(collect(1:length(segs)), box))

    BoxCounter(segs, seg_box, min_segment(curve)/2)
end

# One box intersected. That box is bounding box
start{D, T}(bc::BoxCounter{D, T}) = (1, size(bc))

# Refine and count
function next{D, T}(bc::BoxCounter{D, T}, state)
    
    count = 0
    # Go over all the boxes on the current level
    for i in 1:length(bc.seg_box)

        index_segments, boxes = shift!(bc.seg_box)
        # Are children intersected
        for box in split(boxes)
            isected_segments = Vector{Int}()
            # Collect isects
            for index in index_segments
                collides(bc.segments[index], box) && push!(isected_segments, index)
            end
            # Increate count if there are and send for next round
            if !isempty(isected_segments)
                count += 1
                push!(bc.seg_box, Pair(isected_segments, box))
            end
        end
    end
    # We have halved the mesh
    (count, last(state)/2.)
end

# FIXME: It doesn't make sence to go much beyond the resolution of the curve
done{D, T}(c::BoxCounter{D, T}, state) = true #minimum(last(state)) < c.stop_length


"""Char size of the hyper box - best fit by cuve"""
volume{D, T}(p::SVector{D, T}) = (prod(p))^(1./D)


"""Remember volumes and counts as the refinement is done with box counter"""
mutable struct FractalDim{D, T}
    bc::BoxCounter{D, T}

    volumes::Vector{Float64}
    counts::Vector{Int}
    next_state::Tuple{Int, SVector{D, T}}
end


"""Estimate fd by lin fit of the stepped FractalDim"""
function estimate_fd{D, T}(fd::FractalDim{D, T}, slice=3)
    @assert !isempty(fd.volumes) && !isempty(fd.counts)
    
    if length(fd.counts) < slice
        x = log.(fd.volumes)
        y = log.(fd.counts)
    else
        x = log.(fd.volumes[end-slice:end])
        y = log.(fd.counts[end-slice:end])
    end
    -last(Polynomials.polyfit(x, y, 1).a)
end


"""Curve constructor"""
function FractalDim{D, T}(c::Curve{D, T})
    bc = BoxCounter(c)
    volumes = Vector{Float64}()
    counts = Vector{Int}()
    next_state = (1, size(bc))

    FractalDim(bc, volumes, counts, next_state)
end


"""Start stepping"""
function start{D, T}(fd::FractalDim{D, T})
    state = start(fd.bc)

    push!(fd.counts, first(state))
    push!(fd.volumes, volume(last(state)))
    fd.next_state = state
               
    state
end


"""Evolve"""
function next{D, T}(fd::FractalDim{D, T}, state)
    state = next(fd.bc, state)
    
    push!(fd.counts, first(state))
    push!(fd.volumes, volume(last(state)))
    fd.next_state = state
    
    state
end


"""Let box counter decide"""
done{D, T}(fd::FractalDim{D, T}, state) = done(fd.bc)


"""Evolve the started iterator nlevel times"""
function fractal_dim!{D, T}(fd::FractalDim{D, T}, nlevels::Int)
    state = fd.next_state
    for i in 1:nlevels
        state = next(fd, state)
    end
    (estimate_fd(fd), fd)
end

"""
Step the estimator for curve returning the counter which can be 
further evolved if the precision is not okay
"""
fractal_dim{D, T}(c::Curve{D, T}, nlevels::Int) = fractal_dim!(FractalDim(c), nlevels)
