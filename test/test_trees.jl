# A small capillary stub so the cross-link / leaf-connection searches have targets.
function _stub_network(::Type{T}, n = 6) where {T}
    net = VascularNetwork{T}()
    for k in 1:n
        MM.add_vertex!(net, SVector{3,T}(k, 0, 0))
    end
    return net
end

@testset "add_venule! ($T)" for T in TEST_TYPES
    net = _stub_network(T)
    connecting = Int[]
    MM.add_venule!(net, MersenneTwister(0), NetworkConfig{T}(), one(T), connecting)

    @test eltype(net.vertices) == SVector{3,T}
    @test nvertices(net) == 6 + 11               # 11 venule vertices
    @test nelements(net) == 10 + 5               # axial chain + cross-links
    @test length(connecting) == 5
    @test count(==(Outlet), net.node_tags) == 1
    @test count(==(InletVenule), net.node_tags) == 1
    @test all(>(0), net.radii)
end

@testset "add_arteriole! ($T)" for T in TEST_TYPES
    net = _stub_network(T)
    connecting = Int[]
    MM.add_arteriole!(net, MersenneTwister(0), NetworkConfig{T}(), connecting;
        start = SVector{3,T}(0, 0, 0),
        finish = SVector{3,T}(0, 0, 1),
        trunk_radius = T(9),
        o1 = SVector{3,T}(1, 0, 1),
        o2 = SVector{3,T}(-1, 0, 1),
        o11 = SVector{3,T}(1, 0, 1),
        o12 = SVector{3,T}(2, 0, 1),
        o21 = SVector{3,T}(-1, 0, 1),
        o22 = SVector{3,T}(-2, 0, 1),
        exclude_count = 8)

    @test eltype(net.vertices) == SVector{3,T}
    @test nvertices(net) == 6 + 8                # 8 tree vertices
    @test nelements(net) == 7 + 4               # internal tree + leaf connections
    @test length(connecting) == 4
    @test count(==(Inlet), net.node_tags) == 1
end

@testset "enforce_min_radius! ($T)" for T in TEST_TYPES
    net = VascularNetwork{T}()
    MM.add_vertex!(net, SVector{3,T}(0, 0, 0))
    MM.add_vertex!(net, SVector{3,T}(1, 0, 0))
    MM.add_vertex!(net, SVector{3,T}(2, 0, 0))
    MM.add_element!(net, 1, 2, 0.5)              # incident to vertex 1, thin — widened
    MM.add_element!(net, 1, 3, 99.0)            # incident to vertex 1, already wide — kept
    cfg = NetworkConfig{T}()
    MM.enforce_min_radius!(net, [1], cfg)

    @test net.radii[1] == T(5) * cfg.scale_factor_arteriole
    @test net.radii[2] == 99                    # above the minimum, left untouched
end
