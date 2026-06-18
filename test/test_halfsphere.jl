disk(R) = c -> c[1]^2 + c[2]^2 <= R^2
radius2(v) = v[1]^2 + v[2]^2 + v[3]^2

@testset "dome_frame orthonormal ($T)" for T in TEST_TYPES
    dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
    n, t1, t2 = IM.dome_frame(T, T(0.7), T(1.3))
    @test dot3(n, n) ≈ 1
    @test dot3(t1, t1) ≈ 1
    @test dot3(t2, t2) ≈ 1
    @test abs(dot3(n, t1)) < sqrt(eps(T))
    @test abs(dot3(n, t2)) < sqrt(eps(T))
    @test abs(dot3(t1, t2)) < sqrt(eps(T))
end

@testset "hex_grid_disk_points ($T)" for T in TEST_TYPES
    cfg = NetworkConfig{T}(domain = HalfSphere())
    points, a = IM.hex_grid_disk_points(MersenneTwister(42), cfg)
    @test a isa T
    @test a > 0
    @test !isempty(points)
    @test eltype(points) == NTuple{2,T}
end

@testset "voronoi_network disk predicate ($T)" for T in TEST_TYPES
    cfg = NetworkConfig{T}(domain = HalfSphere())
    R = cfg.sphere_radius
    points, _ = IM.hex_grid_disk_points(MersenneTwister(42), cfg)

    verts, edges = IM.voronoi_network(points, R; in_region = disk(R))
    @test !isempty(verts)
    @test !isempty(edges)
    @test all(v -> v[1]^2 + v[2]^2 <= R^2, verts)

    # the default (square) predicate is unaffected by the new keyword
    sq, _ = IM.voronoi_network(points, R)
    sq2, _ = IM.voronoi_network(points, R)
    @test sq == sq2
end

@testset "extrude_dome! ($T)" for T in TEST_TYPES
    cfg = NetworkConfig{T}(domain = HalfSphere())
    R = cfg.sphere_radius
    points, a = IM.hex_grid_disk_points(MersenneTwister(7), cfg)
    v2d, e2d = IM.voronoi_network(points, R; in_region = disk(R))

    net = VascularNetwork{T}()
    gid, nlayers = IM.extrude_dome!(net, MersenneTwister(0), v2d, e2d, a, R, cfg)

    @test eltype(net.vertices) == SVector{3,T}
    @test eltype(net.radii) == T
    @test size(gid) == (length(v2d), nlayers)
    @test nlayers >= 2

    # layers form a dome: nested active sets shrink from the full flat face to the apex
    active = [count(!=(0), gid[:, ℓ]) for ℓ in 1:nlayers]
    @test active[1] >= length(v2d) - 2                 # flat face keeps (essentially) all
    @test issorted(active; rev = true)
    @test active[1] > active[end]
    @test active[end] >= 1                             # apex non-degenerate

    # gid nonzeros are exactly the global vertex indices, each assigned once
    @test sort(filter(!=(0), vec(gid))) == collect(1:nvertices(net))

    # every element references a present vertex — no dangling off a clipped vertex
    @test all(e -> 1 <= e[1] <= nvertices(net) && 1 <= e[2] <= nvertices(net), net.elements)
    @test all(>(0), net.radii)

    # the bed lies within the dome (allowing for layer jitter) above the flat face
    margin = 6 * max(T(0.1) * a, T(0.005) * R)
    @test all(v -> radius2(v) <= (R + margin)^2, net.vertices)
    @test all(v -> v[3] >= -R / 100, net.vertices)
end

@testset "determinism halfsphere ($T)" for T in TEST_TYPES
    cfg = NetworkConfig{T}(domain = HalfSphere())
    a = generate_network(MersenneTwister(0), cfg)
    b = generate_network(MersenneTwister(0), cfg)
    @test a.vertices == b.vertices
    @test a.elements == b.elements
    @test a.radii == b.radii
    @test a.node_tags == b.node_tags
    @test a.bcvalues == b.bcvalues

    c = generate_network(MersenneTwister(1), cfg)
    @test a.vertices != c.vertices

    @test a isa VascularNetwork{T}
end

@testset "default precision halfsphere" begin
    net = generate_network(MersenneTwister(0), NetworkConfig(domain = HalfSphere()))
    @test net isa VascularNetwork{Float64}
end

@testset "structure halfsphere ($T)" for T in TEST_TYPES
    cfg = NetworkConfig{T}(domain = HalfSphere())
    R = cfg.sphere_radius
    net = generate_network(MersenneTwister(0), cfg)

    @test eltype(net.vertices) == SVector{3,T}
    @test eltype(net.radii) == T
    @test eltype(net.bcvalues) == T

    @test nvertices(net) == length(net.node_tags) == length(net.bcvalues)
    @test nelements(net) == length(net.radii)

    @test all(e -> 1 <= e[1] <= nvertices(net) && 1 <= e[2] <= nvertices(net), net.elements)
    @test all(isfinite, net.radii)
    @test all(>(0), net.radii)
    @test ncomponents(net) == 1                         # no disconnected capillaries

    @test count(==(Inlet), net.node_tags) == 1          # one feeding arteriole
    @test count(==(Outlet), net.node_tags) == 1         # one draining venule
    @test count(==(InletVenule), net.node_tags) == 0    # single-outlet venule
    @test count(==(NoFlow), net.node_tags) > 0          # flat face + dome shell

    @test net.bcvalues[findfirst(==(Outlet), net.node_tags)] == 2000
    @test all(==(T(1e-3)), net.bcvalues[findall(==(Inlet), net.node_tags)])

    # hemisphere-shaped: vessels stay near/inside the dome, above the flat face
    @test all(v -> v[3] >= -R / 100, net.vertices)
    @test all(v -> radius2(v) <= (T(1.2) * R)^2, net.vertices)
end

@testset "config is respected halfsphere ($T)" for T in TEST_TYPES
    sparse = generate_network(MersenneTwister(0), NetworkConfig{T}(domain = HalfSphere(), num_fibers = 40))
    dense = generate_network(MersenneTwister(0), NetworkConfig{T}(domain = HalfSphere(), num_fibers = 90))
    @test nvertices(dense) > nvertices(sparse)

    small = generate_network(MersenneTwister(0), NetworkConfig{T}(domain = HalfSphere(), sphere_radius = 80))
    big = generate_network(MersenneTwister(0), NetworkConfig{T}(domain = HalfSphere(), sphere_radius = 160))
    @test maximum(radius2, big.vertices) > maximum(radius2, small.vertices)
end
