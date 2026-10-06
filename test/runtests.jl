using MicrovascularMeshes
using LinearAlgebra: dot, norm
using Random
using StaticArrays
using Test

const MM = MicrovascularMeshes

# Every type-parametric testset is exercised over each of these element types.
const TEST_TYPES = (Float64, Float32)

# Number of connected components of a network's element graph (1 == fully connected).
function ncomponents(net)
    n = nvertices(net)
    n == 0 && return 0
    parent = collect(1:n)
    function root(x)
        while parent[x] != x
            parent[x] = parent[parent[x]]
            x = parent[x]
        end
        return x
    end
    for e in net.elements
        a, b = root(e[1]), root(e[2])
        a != b && (parent[a] = b)
    end
    return length(Set(root(v) for v in 1:n))
end

@testset "MicrovascularMeshes.jl" begin
    include("test_voronoi.jl")
    include("test_extrude.jl")
    include("test_trees.jl")
    include("test_prune.jl")
    include("test_smooth.jl")
    include("test_toy.jl")
    include("test_generate.jl")
    include("test_halfsphere.jl")
    include("test_scaling.jl")
end
