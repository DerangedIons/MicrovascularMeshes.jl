# A bent 3-vertex path: unit segment along x, then a diagonal to (2, 1, 0).
function _bent_path(::Type{T}) where {T}
    net = VascularNetwork{T}()
    for p in (SVector{3,T}(0, 0, 0), SVector{3,T}(1, 0, 0), SVector{3,T}(2, 1, 0))
        MM.add_vertex!(net, p)
    end
    MM.add_element!(net, 1, 2, one(T))
    MM.add_element!(net, 2, 3, T(2))
    return net
end

@testset "smooth_network invariants and radius inheritance ($T)" for T in TEST_TYPES
    net = _bent_path(T)
    MM.set_bc!(net, 1, Inlet, T(1e-3))
    MM.set_bc!(net, 3, Outlet, T(2000))

    out = smooth_network(MersenneTwister(0), net; max_tilt = T(0), spacing = T(0.25))

    # element 1 has length 1 -> 4 sub-elements; element 2 has length √2 -> 6
    @test nelements(out) == 4 + 6
    @test nvertices(out) == 3 + 3 + 5
    @test out isa VascularNetwork{T}
    @test length(out.radii) == nelements(out)
    @test length(out.node_tags) == nvertices(out)
    @test length(out.bcvalues) == nvertices(out)
    # originals untouched, in place, with their tags and BC values
    @test out.vertices[1:3] == net.vertices
    @test out.node_tags[1:3] == [Inlet, Interior, Outlet]
    @test out.bcvalues[1:3] == [T(1e-3), zero(T), T(2000)]
    # every new vertex is a plain interior sample
    @test all(==(Interior), out.node_tags[4:end])
    @test all(iszero, out.bcvalues[4:end])
    # sub-elements inherit the parent radius
    @test count(==(one(T)), out.radii) == 4
    @test count(==(T(2)), out.radii) == 6
    @test ncomponents(out) == 1
    @test all(e -> 1 <= e[1] <= nvertices(out) && 1 <= e[2] <= nvertices(out), out.elements)
end

@testset "smooth_network keeps a collinear chain exactly straight ($T)" for T in TEST_TYPES
    net = VascularNetwork{T}()
    for p in (SVector{3,T}(0, 0, 0), SVector{3,T}(1, 0, 0), SVector{3,T}(3, 0, 0))
        MM.add_vertex!(net, p)
    end
    MM.add_element!(net, 1, 2, one(T))
    MM.add_element!(net, 2, 3, one(T))

    out = smooth_network(MersenneTwister(0), net; max_tilt = T(0), spacing = T(0.5))

    @test nelements(out) == 2 + 4
    @test nvertices(out) == 3 + 1 + 3
    # all inputs lie on the x-axis, so the Hermite combination is exactly zero off-axis
    @test all(v -> v[2] == 0 && v[3] == 0, out.vertices)
    # chord-scaled tangents reproduce the straight line with uniform parametrization
    @test out.vertices[4] ≈ SVector{3,T}(0.5, 0, 0)
end

@testset "smooth_network is tangent-continuous at a pass-through vertex ($T)" for T in TEST_TYPES
    net = _bent_path(T)
    out = smooth_network(MersenneTwister(0), net; max_tilt = T(0), spacing = T(1 / 16))

    v2 = out.vertices[2]
    into = [e for e in out.elements if e[2] == 2]   # last sub-element of parent 1
    from = [e for e in out.elements if e[1] == 2]   # first sub-element of parent 2
    @test length(into) == 1
    @test length(from) == 1
    u1 = (v2 - out.vertices[only(into)[1]]) / norm(v2 - out.vertices[only(into)[1]])
    u2 = (out.vertices[only(from)[2]] - v2) / norm(out.vertices[only(from)[2]] - v2)
    through = SVector{3,T}(2, 1, 0) / norm(SVector{3,T}(2, 1, 0))  # A -> B direction

    # secants approximate the shared tangent to O(1/n); 0.99 ≈ 8° slack
    @test dot(u1, u2) > 0.99
    @test dot(u1, through) > 0.99
    @test dot(u2, through) > 0.99
end

@testset "smooth_network pairs the most-collinear elements at a junction ($T)" for T in TEST_TYPES
    net = VascularNetwork{T}()
    for p in (SVector{3,T}(-1, 0, 0),   # A
              SVector{3,T}(0, 0, 0),    # C, the junction
              SVector{3,T}(1, 0, 0),    # B
              SVector{3,T}(0, 1, 0))    # D, the branch tip
        MM.add_vertex!(net, p)
    end
    # branch first, so the winning through pair (A-C, C-B) is not first-encountered
    MM.add_element!(net, 2, 4, one(T))
    MM.add_element!(net, 1, 2, one(T))
    MM.add_element!(net, 2, 3, one(T))

    out = smooth_network(MersenneTwister(0), net; max_tilt = T(0), spacing = T(0.25))

    @test out.vertices[1:4] == net.vertices
    # the through vessel stays exactly on the x-axis and the branch exactly on the
    # y-axis, so every sample lies on one of the two axes in the z = 0 plane
    @test all(v -> v[3] == 0, out.vertices)
    @test all(v -> v[1] * v[2] == 0, out.vertices)
    @test ncomponents(out) == 1
end

@testset "smooth_network is deterministic in the rng ($T)" for T in TEST_TYPES
    net = _bent_path(T)
    out1 = smooth_network(MersenneTwister(7), net)
    out2 = smooth_network(MersenneTwister(7), net)
    out3 = smooth_network(MersenneTwister(8), net)

    @test out1.vertices == out2.vertices
    @test out1.elements == out2.elements
    @test out1.radii == out2.radii
    @test out1.node_tags == out2.node_tags
    @test out1.bcvalues == out2.bcvalues
    @test out1.vertices != out3.vertices  # free-end tilts actually depend on the seed
end

@testset "smooth_network passes short elements through unchanged ($T)" for T in TEST_TYPES
    net = _bent_path(T)
    out = smooth_network(MersenneTwister(0), net; spacing = T(100))

    @test out !== net  # always a fresh network, even when nothing is subdivided
    @test out.vertices == net.vertices
    @test out.elements == net.elements
    @test out.radii == net.radii
    @test out.node_tags == net.node_tags
    @test out.bcvalues == net.bcvalues
end

@testset "smooth_network refines a generated network ($T)" for T in TEST_TYPES
    net = generate_network(MersenneTwister(0), NetworkConfig{T}())
    out = smooth_network(MersenneTwister(1), net)

    @test ncomponents(out) == 1
    @test nvertices(out) >= nvertices(net)
    @test nelements(out) >= nelements(net)
    @test length(out.radii) == nelements(out)
    @test length(out.node_tags) == nvertices(out)
    @test length(out.bcvalues) == nvertices(out)
    @test out.vertices[1:nvertices(net)] == net.vertices
    for tag in (Inlet, Outlet, InletVenule, NoFlow)
        @test count(==(tag), out.node_tags) == count(==(tag), net.node_tags)
    end
    @test all(>(0), out.radii)
end

@testset "smooth_network validates arguments" begin
    net = _bent_path(Float64)
    @test_throws ArgumentError smooth_network(MersenneTwister(0), net; spacing = 0.0)
    @test_throws ArgumentError smooth_network(MersenneTwister(0), net; max_tilt = -0.1)
    @test smooth_network(net) isa VascularNetwork{Float64}  # default_rng convenience
end
