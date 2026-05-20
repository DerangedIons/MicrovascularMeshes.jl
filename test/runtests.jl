using IdealizedMicrovasculature
using Random
using StaticArrays
using Test

const IM = IdealizedMicrovasculature

@testset "IdealizedMicrovasculature.jl" begin
    include("test_voronoi.jl")
    include("test_extrude.jl")
    include("test_trees.jl")
    include("test_generate.jl")
end
