@testset "hex_grid_points" begin
    cfg = NetworkConfig()
    points, a = IM.hex_grid_points(MersenneTwister(42), cfg)
    @test a > 0
    @test !isempty(points)
    @test eltype(points) == NTuple{2,Float64}
end

@testset "voronoi_network" begin
    cfg = NetworkConfig()
    points, _ = IM.hex_grid_points(MersenneTwister(42), cfg)
    verts, edges = IM.voronoi_network(points, cfg.region_size)

    @test !isempty(verts)
    @test !isempty(edges)

    # every extracted vertex lies inside the region
    @test all(v -> 0.0 <= v[1] <= cfg.region_size && 0.0 <= v[2] <= cfg.region_size, verts)

    # every edge references valid, distinct vertices
    @test all(e -> 1 <= e[1] <= length(verts) && 1 <= e[2] <= length(verts), edges)
    @test all(e -> e[1] != e[2], edges)

    # tessellating the same points twice gives the same graph
    verts2, edges2 = IM.voronoi_network(points, cfg.region_size)
    @test verts == verts2
    @test edges == edges2
end
