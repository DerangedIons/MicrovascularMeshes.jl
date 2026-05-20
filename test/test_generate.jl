@testset "determinism" begin
    a = generate_network(MersenneTwister(0))
    b = generate_network(MersenneTwister(0))
    @test a.vertices == b.vertices
    @test a.elements == b.elements
    @test a.radii == b.radii
    @test a.node_tags == b.node_tags
    @test a.bcvalues == b.bcvalues

    c = generate_network(MersenneTwister(1))
    @test a.vertices != c.vertices
end

@testset "structure" begin
    net = generate_network(MersenneTwister(0))

    @test nvertices(net) == length(net.node_tags) == length(net.bcvalues)
    @test nelements(net) == length(net.radii)

    @test all(e -> 1 <= e[1] <= nvertices(net) && 1 <= e[2] <= nvertices(net), net.elements)
    @test all(isfinite, net.radii)
    @test all(>(0), net.radii)

    @test count(==(Inlet), net.node_tags) == 2          # two arterioles
    @test count(==(Outlet), net.node_tags) == 1         # one venule
    @test count(==(InletVenule), net.node_tags) == 1
    @test count(==(NoFlow), net.node_tags) > 0          # top + bottom layers

    @test net.bcvalues[findfirst(==(Outlet), net.node_tags)] == 2000.0
    @test net.bcvalues[findfirst(==(InletVenule), net.node_tags)] == 1e-3
    @test all(==(1e-3), net.bcvalues[findall(==(Inlet), net.node_tags)])
end

@testset "config is respected" begin
    sparse = generate_network(MersenneTwister(0), NetworkConfig(num_fibers = 40))
    dense = generate_network(MersenneTwister(0), NetworkConfig(num_fibers = 90))
    @test nvertices(dense) > nvertices(sparse)

    shallow = generate_network(MersenneTwister(0), NetworkConfig(num_layers = 10))
    deep = generate_network(MersenneTwister(0), NetworkConfig(num_layers = 30))
    @test nvertices(deep) > nvertices(shallow)
end
