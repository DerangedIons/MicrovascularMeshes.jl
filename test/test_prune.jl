@testset "prune_disconnected drops orphans ($T)" for T in TEST_TYPES
    net = VascularNetwork{T}()
    for p in (SVector{3,T}(0, 0, 0), SVector{3,T}(1, 0, 0), SVector{3,T}(2, 0, 0),  # main path
              SVector{3,T}(0, 5, 0), SVector{3,T}(1, 5, 0),                          # orphan pair
              SVector{3,T}(0, 9, 0))                                                 # isolated vertex
        MM.add_vertex!(net, p)
    end
    MM.add_element!(net, 1, 2, one(T))
    MM.add_element!(net, 2, 3, one(T))
    MM.add_element!(net, 4, 5, one(T))
    MM.set_bc!(net, 1, Inlet, T(1e-3))

    pruned = MM.prune_disconnected(net)

    @test nvertices(pruned) == 3                       # only the main component survives
    @test nelements(pruned) == 2
    @test ncomponents(pruned) == 1
    @test count(==(Inlet), pruned.node_tags) == 1      # the I/O node is preserved
    @test all(e -> 1 <= e[1] <= 3 && 1 <= e[2] <= 3, pruned.elements)
    @test length(pruned.radii) == nelements(pruned)
    @test length(pruned.bcvalues) == nvertices(pruned)
end

@testset "prune_disconnected keeps connected nets unchanged ($T)" for T in TEST_TYPES
    net = VascularNetwork{T}()
    for p in (SVector{3,T}(0, 0, 0), SVector{3,T}(1, 0, 0), SVector{3,T}(2, 0, 0))
        MM.add_vertex!(net, p)
    end
    MM.add_element!(net, 1, 2, one(T))
    MM.add_element!(net, 2, 3, one(T))

    pruned = MM.prune_disconnected(net)
    @test pruned.vertices == net.vertices
    @test pruned.elements == net.elements
    @test pruned.radii == net.radii
end
