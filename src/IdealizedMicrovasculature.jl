module IdealizedMicrovasculature

using Random: AbstractRNG, default_rng, rand, randn, randperm
using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector
import DelaunayTriangulation as DT

export NodeTag, Interior, NoFlow, Inlet, Outlet, InletVenule
export AbstractDomain, Cube, HalfSphere
export NetworkConfig, VascularNetwork, generate_network, smooth_network, toy_network
export nvertices, nelements
export plotnetwork, plotnetwork!

include("domains.jl")
include("types.jl")
include("voronoi.jl")
include("extrude.jl")
include("trees.jl")
include("prune.jl")
include("smooth.jl")
include("generate.jl")
include("toy.jl")

"""
    plotnetwork(net; kwargs...) -> Figure
    plotnetwork!(ax, net; kwargs...) -> ax

Plot a [`VascularNetwork`](@ref) in 3D: vessels as line segments colored by radius, and
boundary nodes highlighted by [`NodeTag`](@ref).

These functions are implemented by a package extension and require a Makie backend — load
`GLMakie` (interactive) or `CairoMakie` (static) before calling them.

### Keyword Arguments
- `colormap`: colormap for the vessel radii (default `:viridis`).
- `linewidth`: vessel line width (default `1.5`).
- `show_nodes`: whether to scatter the boundary nodes (default `true`).

### Examples
```julia
using GLMakie, IdealizedMicrovasculature
plotnetwork(generate_network())
```
"""
function plotnetwork end
function plotnetwork! end

end # module IdealizedMicrovasculature
