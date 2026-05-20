# A small capillary stub so the cross-link / leaf-connection searches have targets.
function _stub_network(n = 6)
    net = VascularNetwork()
    for k in 1:n
        IM.add_vertex!(net, SVector{3,Float64}(k, 0.0, 0.0))
    end
    return net
end

@testset "add_venule!" begin
    net = _stub_network()
    connecting = Int[]
    IM.add_venule!(net, MersenneTwister(0), NetworkConfig(), 1.0, connecting)

    @test nvertices(net) == 6 + 11               # 11 venule vertices
    @test nelements(net) == 10 + 5               # axial chain + cross-links
    @test length(connecting) == 5
    @test count(==(Outlet), net.node_tags) == 1
    @test count(==(InletVenule), net.node_tags) == 1
    @test all(>(0), net.radii)
end

@testset "add_arteriole!" begin
    net = _stub_network()
    connecting = Int[]
    IM.add_arteriole!(net, MersenneTwister(0), NetworkConfig(), connecting;
        start = SVector{3,Float64}(0.0, 0.0, 0.0),
        finish = SVector{3,Float64}(0.0, 0.0, 1.0),
        trunk_radius = 9.0,
        o1 = SVector{3,Float64}(1.0, 0.0, 1.0),
        o2 = SVector{3,Float64}(-1.0, 0.0, 1.0),
        o11 = SVector{3,Float64}(1.0, 0.0, 1.0),
        o12 = SVector{3,Float64}(2.0, 0.0, 1.0),
        o21 = SVector{3,Float64}(-1.0, 0.0, 1.0),
        o22 = SVector{3,Float64}(-2.0, 0.0, 1.0),
        exclude_count = 8)

    @test nvertices(net) == 6 + 8                # 8 tree vertices
    @test nelements(net) == 7 + 4               # internal tree + leaf connections
    @test length(connecting) == 4
    @test count(==(Inlet), net.node_tags) == 1
end

@testset "enforce_min_radius!" begin
    net = VascularNetwork()
    IM.add_vertex!(net, SVector{3,Float64}(0, 0, 0))
    IM.add_vertex!(net, SVector{3,Float64}(1, 0, 0))
    IM.add_vertex!(net, SVector{3,Float64}(2, 0, 0))
    IM.add_element!(net, 1, 2, 0.5)              # incident to vertex 1, thin — widened
    IM.add_element!(net, 1, 3, 99.0)            # incident to vertex 1, already wide — kept
    cfg = NetworkConfig()
    IM.enforce_min_radius!(net, [1], cfg)

    @test net.radii[1] == 5.0 * cfg.scale_factor_arteriole
    @test net.radii[2] == 99.0                  # above the minimum, left untouched
end
