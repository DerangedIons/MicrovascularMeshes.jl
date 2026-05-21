using IdealizedMicrovasculature
using Random
using StaticArrays
using Test

const IM = IdealizedMicrovasculature

# Every type-parametric testset is exercised over each of these element types.
const TEST_TYPES = (Float64, Float32)

@testset "IdealizedMicrovasculature.jl" begin
    include("test_voronoi.jl")
    include("test_extrude.jl")
    include("test_trees.jl")
    include("test_generate.jl")
end
