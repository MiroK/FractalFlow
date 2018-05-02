using FractalFlow, StaticArrays, PyCall, Polynomials
using FractalFlow: BoundingBox, Segment, BoxCounter

@pyimport matplotlib.pyplot as plt
@pyimport matplotlib.collections as plt_col
@pyimport matplotlib.patches as plt_patch


"""Plot list of lines as line collection"""
function plot{T}(segments::Vector{Segment{2, T}}, ax)
    segments = [[seg.A, seg.B] for seg in segments]
    lines = plt_col.LineCollection(segments, colors="k")
    py"$(ax).add_collection($(lines))"
    ax
end

"""Plot outline of bounding boxes"""
function plot{T}(boxes::Vector{BoundingBox{2, T}}, ax)
    boxes = [plt_patch.Rectangle(box.low, (box.high - box.low)...)
             for box in boxes]
    boxes = plt_col.PatchCollection(boxes, edgecolor="b", facecolor="none")
    py"$(ax).add_collection($(boxes))"
    ax
end


volume{D, T}(p::SVector{D, T}) = (prod(p))^(1./D)


function fd_estimation{T}(c::Curve{2, T}, nlevel::Int)
    counts = []
    sizes = []
    
    bc = BoxCounter(cc)
    state = start(bc)
    @show state

    fig = plt.figure()
    ax = py"$(fig).gca()"

    plot(bc.segments, ax)

    boxes = map(last, bc.seg_box)
    # Extract dim
    min_ = first(boxes).low - 0.1*ones(2)
    max_ = first(boxes).high + 0.1*ones(2)

    xlim = map(first, (min_, max_))
    ylim = map(last, (min_, max_))
    
    plot(boxes, ax)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.axis("off")

    for level in 1:7
        state = next(bc, state)

        fig = plt.figure()
        ax = py"$(fig).gca()"

        plot(bc.segments, ax)
        
        boxes = map(last, bc.seg_box)
        plot(boxes, ax)
        
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.axis("off")
        
        @show state
        push!(counts, first(state))
        push!(sizes, volume(last(state)))
    end

    fig = plt.figure()
    plt.loglog(sizes, counts, "bx-")
    plt.loglog(sizes, sizes.^(-1.)*(last(counts)/(last(sizes)^(-1.))), "r")

    (counts, sizes)
end

# -------------------------------------------------------------------

# segments = [(SVector(0., 0.), SVector(1., 0)),
#             (SVector(1., 0.), SVector(1., -1)),
#             (SVector(1., -1.), SVector(0., -1)),
#             (SVector(0., -1.), SVector(0., 0)),
#             (SVector(0., 0.), SVector(-1., 0)),
#             (SVector(-1., 0.), SVector(-1., 1)),
#             (SVector(-1., 1.), SVector(0., 1)),
#             (SVector(0., 1.), SVector(0., 0)),
#             (SVector(-1., 0), SVector(0., 1)),
#             (SVector(0., -1.), SVector(1., 0))]

# cc = Curve(segments)

cc = koch_q1_curve(6)

(counts, sizes) = fd_estimation(cc, 8)

# Use the last 3 for eastimation
fit = Polynomials.polyfit(log.(counts[end-3:end]), log.(sizes[end-3:end]), 1)
dim = -last(fit.a)
println("Estimated fractal dim $(dim)")

# FIXME: Would be fun to spin this in parallel
dims = [first(FractalFlow.fractal_dim(cc, 10)) for _ in 1:10]
@show dims

dim = mean(dims)
dim_std = std(dims)
@show (dim, dim_std)
