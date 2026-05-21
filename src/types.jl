#--------------------------------------------------------------------------------# Node tags

"""
    NodeTag

Boundary role of a network vertex. Determines whether (and what kind of) a boundary
condition applies to the node.

- `Interior` — interior node, no boundary condition.
- `NoFlow` — no-flow boundary; assigned to the top and bottom extrusion layers.
- `Inlet` — arteriole inflow node; carries a velocity boundary condition.
- `Outlet` — venule outflow node; carries a pressure boundary condition.
- `InletVenule` — venule inflow node; carries a velocity boundary condition.
"""
@enum NodeTag::Int8 begin
    Interior = 0
    NoFlow = 1
    Inlet = 2
    Outlet = 3
    InletVenule = 4
end

#--------------------------------------------------------------------------------# NetworkConfig

"""
    NetworkConfig{T}(; kwargs...)
    NetworkConfig(; kwargs...)

Parameters controlling [`generate_network`](@ref), parametric over the floating-point
element type `T` (default `Float64`). All keyword arguments have defaults matching the
reference idealized network.

Select a non-default precision through the type parameter, e.g.
`NetworkConfig{Float32}(; num_fibers = 60)`. The bare `NetworkConfig(; kwargs...)` form
always yields `NetworkConfig{Float64}` — it cannot infer `T` from the keywords.

### Keyword Arguments
- `num_fibers`: target number of Voronoi seed cells; controls capillary mesh density.
- `region_size`: side length of the square region the capillary mesh covers.
- `random_perturbation`: stddev of the hex-grid jitter, as a fraction of the grid pitch.
- `num_layers`: number of stacked layers the planar mesh is extruded into.
- `layer_spacing_fraction`: spacing between extrusion layers, as a fraction of `region_size`.
- `drop_fraction_hex_grid`: fraction of horizontal (inter-layer) edges dropped per layer.
- `radius_default_small_vessels`: median capillary radius (lognormal location).
- `radius_stddev_small_vessels`: lognormal shape parameter of the capillary radius.
- `scale_factor_venule`: multiplies every venule radius.
- `scale_factor_arteriole`: multiplies every arteriole radius.

### Examples
```julia
cfg = NetworkConfig(num_fibers = 60, num_layers = 25)
cfg32 = NetworkConfig{Float32}(num_fibers = 60)
```
"""
Base.@kwdef struct NetworkConfig{T<:Real}
    num_fibers::Int = 40
    region_size::T = 200.0
    random_perturbation::T = 0.1
    num_layers::Int = 20
    layer_spacing_fraction::T = 0.07
    drop_fraction_hex_grid::T = 0.7
    radius_default_small_vessels::T = 2.6
    radius_stddev_small_vessels::T = 0.15
    scale_factor_venule::T = 1.0
    scale_factor_arteriole::T = 0.7
end

#--------------------------------------------------------------------------------# VascularNetwork

"""
    VascularNetwork{T}

A 3D vascular network as a graph of line elements, parametric over the floating-point
element type `T`.

### Fields
- `vertices`: node positions, one `SVector{3,T}` each.
- `elements`: line elements, each an `SVector{2,Int}` of 1-based vertex indices.
- `radii`: vessel radius per element; `length(radii) == length(elements)`.
- `node_tags`: [`NodeTag`](@ref) per vertex; `length(node_tags) == length(vertices)`.
- `bcvalues`: boundary-condition value per vertex (velocity or pressure, depending on the
  tag); `length(bcvalues) == length(vertices)`. Zero for interior nodes.
"""
struct VascularNetwork{T<:Real}
    vertices::Vector{SVector{3,T}}
    elements::Vector{SVector{2,Int}}
    radii::Vector{T}
    node_tags::Vector{NodeTag}
    bcvalues::Vector{T}
end

VascularNetwork{T}() where {T<:Real} =
    VascularNetwork{T}(SVector{3,T}[], SVector{2,Int}[], T[], NodeTag[], T[])
VascularNetwork() = VascularNetwork{Float64}()

"Number of vertices in the network."
nvertices(net::VascularNetwork) = length(net.vertices)

"Number of line elements in the network."
nelements(net::VascularNetwork) = length(net.elements)

"""
    add_vertex!(net, position) -> Int

Append a vertex at `position` (defaulting to an `Interior` tag and zero BC value) and
return its 1-based index.
"""
function add_vertex!(net::VascularNetwork{T}, position::SVector{3,<:Real}) where {T}
    push!(net.vertices, position)
    push!(net.node_tags, Interior)
    push!(net.bcvalues, zero(T))
    return length(net.vertices)
end

"""
    add_element!(net, i, j, radius) -> Int

Append a line element connecting vertices `i` and `j` with the given `radius` and return
its 1-based index.
"""
function add_element!(net::VascularNetwork, i::Integer, j::Integer, radius::Real)
    push!(net.elements, SVector{2,Int}(i, j))
    push!(net.radii, radius)
    return length(net.elements)
end
