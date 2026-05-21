@testset "drop_horizontal_edges" begin
    rng = MersenneTwister(1)
    @test length(IM.drop_horizontal_edges(rng, 10, 0.0)) == 10
    @test length(IM.drop_horizontal_edges(rng, 10, 1.0)) == 0

    kept = IM.drop_horizontal_edges(rng, 10, 0.5)
    @test length(kept) == 5
    @test allunique(kept)
    @test issorted(kept)
    @test all(i -> 1 <= i <= 10, kept)
end

@testset "extrude_capillary! ($T)" for T in TEST_TYPES
    verts2d = SVector{2,T}[SVector(0, 0), SVector(1, 0), SVector(1, 1), SVector(0, 1)]
    edges2d = [SVector(1, 2), SVector(2, 3), SVector(3, 4)]
    cfg = NetworkConfig{T}(num_layers = 3, drop_fraction_hex_grid = 0.0)

    net = VascularNetwork{T}()
    nvpl = IM.extrude_capillary!(net, MersenneTwister(0), verts2d, edges2d, one(T), cfg)

    @test eltype(net.vertices) == SVector{3,T}
    @test eltype(net.radii) == T
    @test nvpl == 4
    @test nvertices(net) == 3 * nvpl                       # num_layers * nvpl

    # vertical = (num_layers-1)*nvpl ; horizontal = (num_layers-1)*nedges (no drop)
    @test nelements(net) == 2 * nvpl + 2 * length(edges2d)
    @test length(net.radii) == nelements(net)
    @test all(>(0), net.radii)

    # the bottom layer sits exactly at z = 0
    @test all(i -> net.vertices[i][3] == 0, 1:nvpl)
end
