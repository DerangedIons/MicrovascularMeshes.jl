@testset "hex_grid_points ($T)" for T in TEST_TYPES
    cfg = NetworkConfig{T}()
    points, a = MM.hex_grid_points(MersenneTwister(42), cfg)
    @test a isa T
    @test a > 0
    @test !isempty(points)
    @test eltype(points) == NTuple{2,T}
end

@testset "voronoi_network ($T)" for T in TEST_TYPES
    cfg = NetworkConfig{T}()
    points, _ = MM.hex_grid_points(MersenneTwister(42), cfg)
    verts, edges = MM.voronoi_network(points, cfg.region_size)

    @test eltype(verts) == SVector{2,T}
    @test !isempty(verts)
    @test !isempty(edges)

    # every extracted vertex lies inside the region
    @test all(v -> 0 <= v[1] <= cfg.region_size && 0 <= v[2] <= cfg.region_size, verts)

    # every edge references valid, distinct vertices
    @test all(e -> 1 <= e[1] <= length(verts) && 1 <= e[2] <= length(verts), edges)
    @test all(e -> e[1] != e[2], edges)

    # tessellating the same points twice gives the same graph
    verts2, edges2 = MM.voronoi_network(points, cfg.region_size)
    @test verts == verts2
    @test edges == edges2
end
