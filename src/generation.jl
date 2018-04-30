# In general curves that are free of bifurcations can be evolved by inserting
# points between pair. 
function evolve_koch(niters::Int, v0::SVector{2, Float64}, v1::SVector{2, Float64}, evolver::Function)
    @assert v0 != v1
    @assert niters > 0
    # This is now pretty but perhaps faster
    points = Vector{SVector{2, Float64}}([v0, v1])
    while niters > 0
        new_points = Vector{SVector{2, Float64}}()
        for i in 1:(length(points)-1)
            v0, v1 = points[i], points[i+1]
            append!(new_points, [v0, evolver(v0, v1)...])
        end
        push!(new_points, v1)

        points = new_points
        niters -= 1
    end
    points
end


"""Step for the koch curve ---- to  --/\-- """
function koch_t_step(v0::SVector{2, Float64}, v1::SVector{2, Float64}, turn=:right)
    @assert v0 != v1

    v13 = v0*2/3 + v1*1/3
    v23 = v0*1/3 + v1*2/3

    t = (v1 - v0)

    shift = (turn == :left) ? SVector{2, Float64}(t[2], -t[1]) : SVector{2, Float64}(-t[2], t[1])
    
    vmid = 0.5*(v0 + v1) + shift*sqrt(3)/6
    # New points
    [v13, vmid, v23]
end


"""Step for the koch curve ---- to  --∩-- """
function koch_q1_step(v0::SVector{2, Float64}, v1::SVector{2, Float64}, turn=:right)
    @assert v0 != v1

    v13 = v0*2/3 + v1*1/3
    v23 = v0*1/3 + v1*2/3

    t = (v1 - v0)

    shift = (turn == :left) ? SVector{2, Float64}(t[2], -t[1]) : SVector{2, Float64}(-t[2], t[1])
    shift /= 3
    
    # New points
    [v13, v13+shift, v23+shift, v23]
end

"""Step for the koch curve ---- to  --∩u-- """
function koch_q2_step(v0::SVector{2, Float64}, v1::SVector{2, Float64}, turn=:right)
    @assert v0 != v1

    v14 = v0*3/4 + v1*1/4
    v34 = v0*1/4 + v1*3/4
    vmid = v0/2 + v1/2

    t = (v1 - v0)

    shift = (turn == :left) ? SVector{2, Float64}(t[2], -t[1]) : SVector{2, Float64}(-t[2], t[1])
    shift /= 4
    
    # New points
    [v14, v14+shift, vmid+shift, vmid, vmid-shift, v34-shift, v34]
end



for evolve in (:koch_t_step, :koch_q1_step, :koch_q2_step)
    curve_name = Symbol(join([split(string(evolve), "_")[1:end-1]..., "curve"], "_"))

    @eval begin
        function $(curve_name)(niters::Int,
                               v0::SVector{2, Float64}=SVector(0., 0.),
                               v1::SVector{2, Float64}=SVector(1., 0.))
            Curve(evolve_koch(niters, v0, v1, $(evolve)))
        end
    end
end


# TODO: hilbert curve (space filling)
#       dragon curve
