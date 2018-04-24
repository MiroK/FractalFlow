using StaticArrays

@test begin
    x0 = SVector(0., 0.)
    x1 = SVector(1., 1.)

    curve = Curve([x0, x1])
    is_consistent(curve)
end

@test begin
    x0 = SVector(0., 0.)
    x1 = SVector(1., 1.)
    x2 = SVector(1., 1.)

    curve = Curve([x0, x1, x2])
    ! is_consistent(curve)
end

@test begin
    x0 = SVector(0., 0.)
    x1 = SVector(0., 0.)

    curve = Curve([x0, x1])
    !is_consistent(curve)
end

@test begin
    x0 = SVector(0., 0.)
    x1 = SVector(1., 1.)
    x2 = SVector(2., 2.)

    curve = Curve([x0, x1, x2], [(1, 2), (2, 3), (3, 1)])
    is_consistent(curve)
end

@test begin
    x0 = SVector(0., 0.)
    x1 = SVector(1., 1.)
    x2 = SVector(2., 2.)

    curve = Curve([x0, x1, x2], [(1, 1), (2, 3), (3, 1)])
    !is_consistent(curve)
end

@test begin
    x0 = SVector(0., 0.)
    x1 = SVector(1., 1.)
    x2 = SVector(2., 2.)

    curve = Curve([x0, x1, x2], [(1, 9), (2, 3), (3, 1)])
    !is_consistent(curve)
end

# IO
@test begin
    ext = ".vtu"
    all(length(Curve(vtu_path*ext).points) > 0 for vtu_path in ("mesh2d", "mesh2din3d", "tet_mesh"))
end

@test begin
    ext = ".xml"
    all(length(Curve(vtu_path*ext).points) > 0 for vtu_path in ("mesh2d", "mesh2din3d", "tet_mesh"))
end

@test begin
    ext = ".msh"
    all(length(Curve(vtu_path*ext).points) > 0 for vtu_path in ("mesh2d", "mesh2din3d", "tet_mesh"))
end
