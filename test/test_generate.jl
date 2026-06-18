@testset "determinism ($T)" for T in TEST_TYPES
    a = generate_network(MersenneTwister(0), NetworkConfig{T}())
    b = generate_network(MersenneTwister(0), NetworkConfig{T}())
    @test a.vertices == b.vertices
    @test a.elements == b.elements
    @test a.radii == b.radii
    @test a.node_tags == b.node_tags
    @test a.bcvalues == b.bcvalues

    c = generate_network(MersenneTwister(1), NetworkConfig{T}())
    @test a.vertices != c.vertices

    @test a isa VascularNetwork{T}
end

@testset "default precision" begin
    @test generate_network(MersenneTwister(0)) isa VascularNetwork{Float64}
end

@testset "structure ($T)" for T in TEST_TYPES
    net = generate_network(MersenneTwister(0), NetworkConfig{T}())

    @test eltype(net.vertices) == SVector{3,T}
    @test eltype(net.radii) == T
    @test eltype(net.bcvalues) == T

    @test nvertices(net) == length(net.node_tags) == length(net.bcvalues)
    @test nelements(net) == length(net.radii)

    @test all(e -> 1 <= e[1] <= nvertices(net) && 1 <= e[2] <= nvertices(net), net.elements)
    @test all(isfinite, net.radii)
    @test all(>(0), net.radii)
    @test ncomponents(net) == 1                         # no disconnected capillaries

    @test count(==(Inlet), net.node_tags) == 2          # two arterioles
    @test count(==(Outlet), net.node_tags) == 1         # one venule
    @test count(==(InletVenule), net.node_tags) == 1
    @test count(==(NoFlow), net.node_tags) > 0          # top + bottom layers

    @test net.bcvalues[findfirst(==(Outlet), net.node_tags)] == 2000
    @test net.bcvalues[findfirst(==(InletVenule), net.node_tags)] == T(1e-3)
    @test all(==(T(1e-3)), net.bcvalues[findall(==(Inlet), net.node_tags)])
end

@testset "config is respected ($T)" for T in TEST_TYPES
    sparse = generate_network(MersenneTwister(0), NetworkConfig{T}(num_fibers = 40))
    dense = generate_network(MersenneTwister(0), NetworkConfig{T}(num_fibers = 90))
    @test nvertices(dense) > nvertices(sparse)

    shallow = generate_network(MersenneTwister(0), NetworkConfig{T}(num_layers = 10))
    deep = generate_network(MersenneTwister(0), NetworkConfig{T}(num_layers = 30))
    @test nvertices(deep) > nvertices(shallow)
end
