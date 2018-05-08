mutable struct Fractal{D, T, N, M}
    segments::Deque{Segment{D, T}}
    step_in::MVector{N, Segment{D, T}}
    step_out::MVector{M, Segment{D, T}}
    step!::Function

    state::Int
end


function Fractal{D, T}(input::Vector{SVector{D, T}}, M::Int, step!::Function)
    @assert length(input) > 1
    
    segments = Deque{Segment{D, T}}()
    for (v0, v1) in zip(input[1:end-1], input[2:end])
        push!(segments, Segment(v0, v1))
    end

    N = length(input) - 1
    step_in = MVector((Segment(SVector(zeros(T, D)...), SVector(zeros(T, D)...)) for i in 1:N)...)

    step_out = MVector((Segment(SVector(zeros(T, D)...), SVector(zeros(T, D)...)) for i in 1:M)...)

    Fractal(segments, step_in, step_out, step!, length(segments))
end

start{D, T, N, M}(f::Fractal{D, T, N, M}) = f.state

function next{D, T, N, M}(f::Fractal{D, T, N, M}, state::Int)
    
    taken, state = state, 0
    while taken > 0
        taken -= N

        for i in 1:N
            f.step_in[i] = shift!(f.segments)
        end
        f.step!(f.step_in, f.step_out)

        for j in 1:M
            push!(f.segments, f.step_out[j])
        end
        state += M
    end
    @assert length(f.segments) == state
    f.state = state
    state
end


done{D, T, N, M}(f::Fractal{D, T, N, M}) = true


function refine!{D, T, N, M}(f::Fractal{D, T, N, M}, count::Int)
    @assert count >= 1
    next(f, f.state)
    (count == 1) ? f : refine!(f, count-1)
end


function koch_t!{T}(input::MVector{1, Segment{2, T}}, output::MVector{4, Segment{2, T}})
    v0, v1 = input[1].A, input[1].B

    v13 = v0*2/3 + v1*1/3
    v23 = v0*1/3 + v1*2/3

    t = (v1 - v0)

    shift = SVector(t[2], -t[1])
    
    vmid = 0.5*(v0 + v1) + shift*sqrt(3)/6
    
    # New points
    output[1] = Segment(v0, v13)
    output[2] = Segment(v13, vmid)
    output[3] = Segment(vmid, v23)
    output[4] = Segment(v23, v1)
end

function KochTFractal{T}(A::SVector{2, T}=SVector(0., 0.), B::SVector{2, T}=SVector(1., 0.))
    Fractal([A, B], 4, koch_t!)
end

# --------------------------------------------------------------------

function koch_q1!{T}(input::MVector{1, Segment{2, T}}, output::MVector{5, Segment{2, T}})
    v0, v1 = input[1].A, input[1].B

    v13 = v0*2/3 + v1*1/3
    v23 = v0*1/3 + v1*2/3

    t = (v1 - v0)
    shift = SVector{2, Float64}(t[2], -t[1])
    shift /= 3
    
    # New points
    #[v13, v13+shift, v23+shift, v23]
    shift_v13 = v13 + shift
    shift_v23 = v23 + shift

    output[1] = Segment(v0, v13)
    output[2] = Segment(v13, shift_v13)
    output[3] = Segment(shift_v13, shift_v23)
    output[4] = Segment(shift_v23, v23)
    output[5] = Segment(v23, v1)
end

function KochQ1Fractal{T}(A::SVector{2, T}=SVector(0., 0.), B::SVector{2, T}=SVector(1., 0.))
    Fractal([A, B], 5, koch_q1!)
end

# --------------------------------------------------------------------

function koch_q2!{T}(input::MVector{1, Segment{2, T}}, output::MVector{8, Segment{2, T}})
    v0, v1 = input[1].A, input[1].B

    v14 = v0*3/4 + v1*1/4
    v34 = v0*1/4 + v1*3/4
    vmid = v0/2 + v1/2

    t = (v1 - v0)

    shift = SVector{2, Float64}(t[2], -t[1])
    shift /= 4
    
    shift_v14 = v14 + shift
    up_vmid = vmid + shift
    down_vmid = vmid - shift
    shift_v34 = v34 - shift
    
    # New points
    # [v14, v14+shift, vmid+shift, vmid, vmid-shift, v34-shift, v34]
    output[1] = Segment(v0, v14)
    output[2] = Segment(v14, shift_v14)
    output[3] = Segment(shift_v14, up_vmid)
    output[4] = Segment(up_vmid, vmid)
    output[5] = Segment(vmid, down_vmid)
    output[6] = Segment(down_vmid, shift_v34)
    output[7] = Segment(shift_v34, v34)
    output[8] = Segment(v34, down_v1)

end

function KochQ2Fractal{T}(A::SVector{2, T}=SVector(0., 0.), B::SVector{2, T}=SVector(1., 0.))
    Fractal([A, B], 8, koch_q2!)
end

# --------------------------------------------------------------------

function levy!{T}(input::MVector{1, Segment{2, T}}, output::MVector{2, Segment{2, T}})
    v0, v1 = input[1].A, input[1].B

    t = (v1 - v0)

    shift = SVector(t[2], -t[1])
    
    vmid = 0.5*(v0 + v1) + shift/2
    
    # New points
    output[1] = Segment(v0, vmid)
    output[2] = Segment(vmid, v1)
end


function LevyFractal{T}(A::SVector{2, T}=SVector(0., 0.), B::SVector{2, T}=SVector(1., 0.))
    Fractal([A, B], 2, levy!)
end

# --------------------------------------------------------------------

function dragon!{T}(input::MVector{2, Segment{2, T}}, output::MVector{4, Segment{2, T}})
    v0, v1, v2 = input[1].A, input[1].B, input[2].B

    v01 = v0/2 + v1/2
    v12 = v1/2 + v2/2

    t1 = (v1 - v0)
    v01 += SVector(t1[2], -t1[1])/2

    t2 = (v2 - v1)
    v12 -= SVector(t2[2], -t2[1])/2

    output[1] = Segment(v0, v01)
    output[2] = Segment(v01, v1)
    output[3] = Segment(v1, v12)
    output[4] = Segment(v12, v2)
end

function DragonFractal{T}(A::SVector{2, T}=SVector(0., 0.),
                          B::SVector{2, T}=SVector(0.5, -0.5),
                          C::SVector{2, T}=SVector(1., 0.))
    Fractal([A, B, C], 4, dragon!)
end

# --------------------------------------------------------------------

function sierpinski_arrow_head!{T}(input::MVector{3, Segment{2, T}}, output::MVector{9, Segment{2, T}})
    v0, v1, v2, v3 = input[1].A, input[2].A, input[3].A, input[3].B

    vmid = (v0 + v3)/2
    
    t = (v3 - v0)

    v01 = (v0 + vmid)/2
    v10 = (v1 + vmid)/2

    v23 = (v2 + vmid)/2
    v32 = (v3 + vmid)/2

    vtop = vmid + 2*((v1+v2)/2 - vmid)
    
    v12 = (v1 + vtop)/2
    v21 = (v2 + vtop)/2

    output[1] = Segment(v0, v01)
    output[2] = Segment(v01, v10)
    output[3] = Segment(v10, v1)
    output[4] = Segment(v1, v12)
    output[5] = Segment(v12, v21)
    output[6] = Segment(v21, v2)
    output[7] = Segment(v2, v23)
    output[8] = Segment(v23, v32)
    output[9] = Segment(v32, v3)
end


function SierpinskiArrowHeadFractal{T}(A::SVector{2, T}=SVector(0., 0.),
                                       B::SVector{2, T}=SVector(0.25, sqrt(3)/4.),
                                       C::SVector{2, T}=SVector(0.75, sqrt(3)/4.),
                                       D::SVector{2, T}=SVector(1., 0.))
    Fractal([A, B, C, D], 9, sierpinski_arrow_head!)
end

# --------------------------------------------------------------------

function fractal_curve{D, T, N, M}(f::Fractal{D, T, N, M})
    segments = Vector{Segment{D, T}}(length(f.segments))
    for (i, seg) in enumerate(f.segments)
        segments[i] = seg
    end

    Curve(segments)
end

Curve{D, T, N, M}(f::Fractal{D, T, N, M}) = Curve(collect(f.segments))
