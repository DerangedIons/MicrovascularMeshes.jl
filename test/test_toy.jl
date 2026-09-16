const TOY_LAYOUTS = (:volumetric, :flat, :radial)

# Closed-form sizes from the toy_network docstring, with n = 2^G leaves per tree;
# `ncap` is the number of capillary elements.
function _toy_expected(layout, G; copies=4)
    n = 2^G
    layout ≡ :radial || return (; nv=4n, ne=6n - 3, ncap=2n - 1)
    side_v = 2 + copies * (2n - 2)
    side_e = 1 + copies * (2n - 2)
    return (; nv=2side_v, ne=2side_e + copies * (2n - 1), ncap=copies * (2n - 1))
end

# Element indices incident to each vertex.
function _incident(net)
    inc = [Int[] for _ in 1:nvertices(net)]
    for (k, e) in enumerate(net.elements)
        push!(inc[e[1]], k)
        push!(inc[e[2]], k)
    end
    return inc
end

# Indices of elements of `net` that `f` does NOT map onto some element of `net` with the
# same radius (undirected). Empty means `f` is a symmetry of the element set.
function _asymmetric_elements(net, f; atol)
    segs = map(k -> (net.vertices[net.elements[k][1]], net.vertices[net.elements[k][2]],
        net.radii[k]), eachindex(net.elements))
    matches((p, q, r), (p′, q′, r′)) = isapprox(r, r′) &&
        ((isapprox(p, p′; atol) && isapprox(q, q′; atol)) ||
         (isapprox(p, q′; atol) && isapprox(q, p′; atol)))
    return filter(eachindex(segs)) do k
        (p, q, r) = segs[k]
        !any(s -> matches((f(p), f(q), r), s), segs)
    end
end

# Vertices of `net` that `f` does NOT map onto some vertex of `net`.
_asymmetric_vertices(net, f; atol) =
    filter(v -> !any(w -> isapprox(f(v), w; atol), net.vertices), net.vertices)

# Rotation by θ about the x axis.
_rotate_x(θ) = p -> typeof(p)(p[1], p[2] * cos(θ) - p[3] * sin(θ), p[2] * sin(θ) + p[3] * cos(θ))

# The half-turn about the bed center that maps every toy layout onto itself: arterial leaf
# k ↦ venous leaf n+1-k, which carries the chiral ladder onto itself.
_half_turn(extent) = p -> typeof(p)(extent - p[1], -p[2], -p[3])

#--------------------------------------------------------------------------------# Structure

@testset "toy_network sizes and graph invariants (layout = $layout, T = $T)" for layout in TOY_LAYOUTS, T in TEST_TYPES
    for G in 1:3
        @testset "generations = $G" begin
            net = toy_network(T; layout, generations=G)
            want = _toy_expected(layout, G)
            (nv, ne) = (nvertices(net), nelements(net))

            @test net isa VascularNetwork{T}
            @test eltype(net.vertices) == SVector{3,T}
            @test nv == want.nv
            @test ne == want.ne
            @test length(net.radii) == ne
            @test length(net.node_tags) == nv
            @test length(net.bcvalues) == nv

            # capillaries are exactly the elements at capillary_radius; tree vessels are wider
            @test count(==(T(2.6)), net.radii) == want.ncap
            @test all(≥(T(2.6)), net.radii)

            # well-formed graph: indices in range, no self-loops, no duplicate elements,
            # no two vertices at the same position, one connected component
            @test all(e -> all(i -> 1 ≤ i ≤ nv, e), net.elements)
            @test isempty(filter(e -> e[1] == e[2], net.elements))
            @test allunique(map(e -> minmax(e[1], e[2]), net.elements))
            @test allunique(net.vertices)
            @test ncomponents(net) == 1
        end
    end
end

@testset "toy_network tags and boundary values (layout = $layout, T = $T)" for layout in TOY_LAYOUTS, T in TEST_TYPES
    extent = 200
    net = toy_network(T; layout, extent)
    deg = length.(_incident(net))
    inlets = findall(==(Inlet), net.node_tags)
    outlets = findall(==(Outlet), net.node_tags)
    others = findall(==(Interior), net.node_tags)

    # one Inlet/Outlet pair, nothing else tagged
    @test length(inlets) == 1
    @test length(outlets) == 1
    @test length(others) == nvertices(net) - 2
    @test all(==(T(1e-3)), net.bcvalues[inlets])
    @test all(==(T(2000)), net.bcvalues[outlets])
    @test all(iszero, net.bcvalues[others])

    # roots are degree-1 trunk stubs on the ends of the x range
    @test all(==(1), deg[inlets])
    @test all(==(1), deg[outlets])
    @test all(v -> net.vertices[v][1] == 0, inlets)
    @test all(v -> net.vertices[v][1] == extent, outlets)
    @test all(v -> 0 < net.vertices[v][1] < extent, others)
end

@testset "toy_network obeys Murray's law at tree junctions (layout = $layout, T = $T)" for layout in TOY_LAYOUTS, T in TEST_TYPES
    (G, copies) = (3, 4)
    rc = T(2.6)
    net = toy_network(T; layout, generations=G, copies)
    inc = _incident(net)
    is_capillary(k) = net.radii[k] == rc
    # r_parent³ - Σ r_child³, taking the widest incident vessel as the parent
    function murray_residual(v)
        rs = net.radii[inc[v]]
        return 2 * maximum(rs)^3 - sum(r -> r^3, rs)
    end

    # tree junctions: branching vertices that touch no capillary (hubs and internal nodes)
    junctions = filter(v -> length(inc[v]) ≥ 3 && !any(is_capillary, inc[v]), 1:nvertices(net))
    fans = layout ≡ :radial ? copies : 1
    @test length(junctions) == 2 * (1 + fans * (2^G - 2))
    bad = filter(v -> !isapprox(murray_residual(v), 0; atol=sqrt(eps(T)) * maximum(net.radii)^3), junctions)
    @test isempty(bad)

    # trunk carries r₀ = rc·2^((G+1)/3), scaled by cbrt(copies) where the fans share a hub
    trunks = map(v -> only(inc[v]), findall(!=(Interior), net.node_tags))
    @test all(k -> net.radii[k] ≈ rc * T(2)^(T(G + 1) / 3) * cbrt(T(fans)), trunks)

    # each ladder has two single-capillary end leaves; the other leaves split exactly into
    # two capillaries
    leaves = filter(v -> any(is_capillary, inc[v]), 1:nvertices(net))
    @test length(leaves) == fans * 2 * 2^G
    @test count(v -> length(inc[v]) == 2, leaves) == 2 * fans
    full = filter(v -> length(inc[v]) == 3, leaves)
    @test length(full) == fans * 2 * (2^G - 1)
    @test isempty(filter(v -> !isapprox(murray_residual(v), 0; atol=sqrt(eps(T)) * rc^3), full))
end

#--------------------------------------------------------------------------------# Geometry

@testset "toy_network flat geometry ($T)" for T in TEST_TYPES
    (G, extent, width) = (3, T(200), T(100))
    n = 2^G
    net = toy_network(T; layout=:flat, generations=G, extent, width)
    inc = _incident(net)
    leaves = filter(v -> any(k -> net.radii[k] == T(2.6), inc[v]), 1:nvertices(net))
    arterial = sort(map(v -> net.vertices[v][2], filter(v -> net.vertices[v][1] < extent / 2, leaves)))
    venous = sort(map(v -> net.vertices[v][2], filter(v -> net.vertices[v][1] > extent / 2, leaves)))

    @test all(v -> v[3] == 0, net.vertices)
    @test all(v -> abs(v[2]) < width / 2, net.vertices)
    # leaves sit on the bed edges, evenly spaced at pitch width/n and centered on y = 0
    @test all(v -> net.vertices[v][1] ≈ T(0.4) * extent || net.vertices[v][1] ≈ T(0.6) * extent, leaves)
    @test length(arterial) == n
    @test arterial ≈ collect(range(-(width / 2) * (1 - T(1) / n), (width / 2) * (1 - T(1) / n); length=n))
    @test venous ≈ arterial
    # hubs at 0.15 and 0.85 of the extent
    hubs = filter(v -> net.vertices[v][2] == 0 && net.node_tags[v] == Interior, 1:nvertices(net))
    @test sort(map(v -> net.vertices[v][1], hubs)) ≈ [T(0.15) * extent, T(0.85) * extent]
end

@testset "toy_network volumetric geometry (G = $G, $T)" for G in 1:4, T in TEST_TYPES
    width = T(100)
    net = toy_network(T; layout=:volumetric, generations=G, width)
    inc = _incident(net)
    leaves = filter(v -> any(k -> net.radii[k] == T(2.6), inc[v]) && net.vertices[v][1] < 100,
        1:nvertices(net))
    ys = unique(map(v -> net.vertices[v][2], leaves))
    zs = unique(map(v -> net.vertices[v][3], leaves))

    # leaves fill a 2^⌈G/2⌉ × 2^⌊G/2⌋ grid in y–z, all inside the width × width square
    @test length(leaves) == 2^G
    @test length(ys) == 2^cld(G, 2)
    @test length(zs) == 2^fld(G, 2)
    @test allunique(map(v -> (net.vertices[v][2], net.vertices[v][3]), leaves))
    @test all(v -> abs(v[2]) < width / 2 && abs(v[3]) < width / 2, net.vertices)
    if G == 1
        # one split means no z split: identical to :flat
        flat = toy_network(T; layout=:flat, generations=G)
        @test net.vertices == flat.vertices
        @test net.elements == flat.elements
        @test net.radii == flat.radii
    else
        @test length(zs) > 1
    end
end

@testset "toy_network radial geometry (copies = $copies, $T)" for copies in 1:5, T in TEST_TYPES
    G = 2
    net = toy_network(T; layout=:radial, generations=G, copies)
    inc = _incident(net)
    axis = filter(v -> v[2] == 0 && v[3] == 0, net.vertices)
    off_axis = filter(v -> !(v[2] == 0 && v[3] == 0), net.vertices)
    # azimuth folded into [0, π): each fan lies in one plane through the x axis
    planes = map(v -> mod(atan(v[3], v[2]), T(π)), off_axis)
    want_planes = map(c -> T(π) * (c - 1) / copies, 1:copies)
    hub_degrees = sort(map(length, inc[findall(v -> v[2] == 0 && v[3] == 0 && 0 < v[1] < 200, net.vertices)]))

    @test length(axis) == 4  # inlet, outlet, and the two hubs
    @test hub_degrees == [2copies + 1, 2copies + 1]
    @test all(φ -> any(φ′ -> isapprox(φ, φ′; atol=1e-4) || isapprox(φ, φ′ + T(π); atol=1e-4), want_planes), planes)
    @test all(φ′ -> count(φ -> isapprox(φ, φ′; atol=1e-4) || isapprox(φ, φ′ + T(π); atol=1e-4), planes) == 2 * (2^(G + 1) - 2), want_planes)

    # rotating by π/copies carries the vertex set and the trees onto themselves; the ladder
    # is chiral, so the capillaries are not asserted here (see the half-turn testset)
    rot = _rotate_x(T(π) / copies)
    @test isempty(_asymmetric_vertices(net, rot; atol=1e-3))
    trees = VascularNetwork{T}(net.vertices, net.elements[net.radii .> T(2.6)],
        net.radii[net.radii .> T(2.6)], net.node_tags, net.bcvalues)
    @test isempty(_asymmetric_elements(trees, rot; atol=1e-3))

    if copies == 1
        flat = toy_network(T; layout=:flat, generations=G)
        @test net.vertices == flat.vertices
        @test net.elements == flat.elements
        @test net.radii == flat.radii
        @test net.node_tags == flat.node_tags
    end
end

@testset "toy_network symmetries (layout = $layout, T = $T)" for layout in TOY_LAYOUTS, T in TEST_TYPES
    extent = T(200)
    for G in 1:3
        net = toy_network(T; layout, generations=G, extent)
        mirror = p -> typeof(p)(extent - p[1], p[2], p[3])
        # the venous tree is the x-mirror of the arterial tree
        @test isempty(_asymmetric_vertices(net, mirror; atol=1e-3))
        # with the chiral ladder, the whole element set is symmetric under a half-turn
        @test isempty(_asymmetric_elements(net, _half_turn(extent); atol=1e-3))
    end
    # and the ladder really is chiral: a pure x-mirror does not preserve the elements
    @test !isempty(_asymmetric_elements(toy_network(T; layout, extent), p -> typeof(p)(extent - p[1], p[2], p[3]); atol=1e-3))
end

@testset "toy_network scales with extent and width ($T)" for T in TEST_TYPES
    net = toy_network(T; layout=:flat, extent=10, width=4, capillary_radius=0.5, generations=2)
    xs = map(v -> v[1], net.vertices)
    ys = map(v -> v[2], net.vertices)

    @test extrema(xs) == (0, 10)
    @test maximum(abs, ys) ≈ T(4) / 2 * (1 - T(1) / 4)
    @test count(==(T(0.5)), net.radii) == 2 * 4 - 1
    # radii are linear in capillary_radius
    @test net.radii ≈ toy_network(T; layout=:flat, generations=2).radii .* (T(0.5) / T(2.6))
end

@testset "toy_network defaults to the volumetric layout ($T)" for T in TEST_TYPES
    net = toy_network(T)
    volumetric = toy_network(T; layout=:volumetric)

    @test net.vertices == volumetric.vertices
    @test net.elements == volumetric.elements
    @test net.radii == volumetric.radii
    @test nvertices(net) == 32
    @test nelements(net) == 45
    @test any(v -> v[3] ≠ 0, net.vertices)  # genuinely 3D, not the flat picture
end

@testset "toy_network is deterministic and precision-preserving (layout = $layout)" for layout in TOY_LAYOUTS
    a = toy_network(; layout)
    b = toy_network(; layout)
    f32 = toy_network(Float32; layout)

    @test a isa VascularNetwork{Float64}
    @test a.vertices == b.vertices
    @test a.elements == b.elements
    @test a.radii == b.radii
    @test f32 isa VascularNetwork{Float32}
    @test eltype(f32.radii) == Float32
    @test eltype(f32.bcvalues) == Float32
    @test f32.elements == a.elements
    @test f32.node_tags == a.node_tags
    @test map(v -> Float32.(v), a.vertices) ≈ f32.vertices
end

#--------------------------------------------------------------------------------# Composition

@testset "toy_network composes with smooth_network (layout = $layout, T = $T)" for layout in TOY_LAYOUTS, T in TEST_TYPES
    net = toy_network(T; layout)
    out = smooth_network(MersenneTwister(0), net)
    nv = nvertices(net)

    @test out isa VascularNetwork{T}
    @test ncomponents(out) == 1
    @test nelements(out) > nelements(net)
    @test out.vertices[1:nv] == net.vertices
    @test out.node_tags[1:nv] == net.node_tags
    @test out.bcvalues[1:nv] == net.bcvalues
    for tag in (Inlet, Outlet)
        @test count(==(tag), out.node_tags) == count(==(tag), net.node_tags)
    end
    # junction degrees are preserved: subdivision only adds degree-2 samples
    @test length.(_incident(out))[1:nv] == length.(_incident(net))
    @test all(==(2), length.(_incident(out))[nv+1:end])
end

#--------------------------------------------------------------------------------# Validation

@testset "toy_network validates arguments" begin
    @test_throws ArgumentError toy_network(; layout=:bogus)
    @test_throws ArgumentError toy_network(; layout="flat")
    @test_throws ArgumentError toy_network(; generations=0)
    @test_throws ArgumentError toy_network(; layout=:stacked)  # removed layout
    @test_throws ArgumentError toy_network(; layout=:radial, copies=0)
    @test_throws ArgumentError toy_network(; extent=0)
    @test_throws ArgumentError toy_network(; width=-1)
    @test_throws ArgumentError toy_network(; capillary_radius=0)
    @test_throws ArgumentError toy_network(Float32; extent=-1.0)
end
