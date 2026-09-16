#--------------------------------------------------------------------------------# Toy network

const _TOY_LAYOUTS = (:volumetric, :flat, :radial)

# Split axis of every generation, one vector per fan the layout hangs off a hub. `:radial`
# spaces its planes by π/copies, not 2π/copies: a fan is symmetric under y ↦ -y, so a
# half-turn maps it onto itself and an even `copies` would stack fans on top of each other.
function _fan_axes(layout::Symbol, ::Type{T}, generations::Integer, copies::Integer) where {T}
    (ŷ, ẑ) = (SVector{3,T}(0, 1, 0), SVector{3,T}(0, 0, 1))
    layout ≡ :volumetric && return [map(g -> isodd(g) ? ŷ : ẑ, 1:generations)]
    layout ≡ :radial || return [fill(ŷ, generations)]
    return map(1:copies) do c
        φ = T(π) * (c - 1) / copies
        fill(SVector{3,T}(0, cos(φ), sin(φ)), generations)
    end
end

# Transverse offsets of every node in one bifurcating fan, grouped by generation:
# `levels[g+1][j]` is node j of generation g (the hub is generation 0), and its children
# are nodes 2j-1 and 2j, displaced by ∓δ along `axes[g]`. δ starts at width/4 and halves
# each time an axis is reused, so every node sits at the centroid of its descendant leaves.
function _fan_levels(axes::Vector{SVector{3,T}}, width::T) where {T}
    levels = [[zero(SVector{3,T})]]
    for (g, â) in enumerate(axes)
        δ = width / 2^(count(==(â), view(axes, 1:g)) + 1)
        push!(levels, mapreduce(p -> [p - δ * â, p + δ * â], vcat, last(levels)))
    end
    return levels
end

# Append the nodes of one fan below `hub` and return its leaf indices. `xs[g+1]` is the x
# of generation g and `radii[g+1]` the radius of the vessel entering a generation-g node.
function _add_fan!(
    net::VascularNetwork{T},
    hub::Int,
    levels::Vector{Vector{SVector{3,T}}},
    xs::Vector{T},
    radii::Vector{T},
) where {T}
    parents = [hub]
    for g in 1:length(levels)-1
        children = Vector{Int}(undef, length(levels[g+1]))
        for (j, offset) in enumerate(levels[g+1])
            children[j] = add_vertex!(net, SVector{3,T}(xs[g+1], 0, 0) + offset)
            add_element!(net, parents[(j+1)>>1], children[j], radii[g+1])
        end
        parents = children
    end
    return parents
end

# Append a root → hub trunk plus every fan below the shared hub; return the root index and
# each fan's leaf indices.
function _add_tree!(
    net::VascularNetwork{T},
    fans::Vector{Vector{Vector{SVector{3,T}}}},
    xroot::T,
    xs::Vector{T},
    radii::Vector{T},
    trunk_radius::T,
) where {T}
    root = add_vertex!(net, SVector{3,T}(xroot, 0, 0))
    hub = add_vertex!(net, SVector{3,T}(xs[1], 0, 0))
    add_element!(net, root, hub, trunk_radius)
    leaves = map(levels -> _add_fan!(net, hub, levels, xs, radii), fans)
    return (root, leaves)
end

# Anastomotic capillary ladder: arterial leaf k feeds venous leaves k and k+1.
function _add_ladder!(net::VascularNetwork, arterial::Vector{Int}, venous::Vector{Int}, radius)
    for k in eachindex(arterial)
        add_element!(net, arterial[k], venous[k], radius)
        k < lastindex(venous) && add_element!(net, arterial[k], venous[k+1], radius)
    end
    return nothing
end

"""
    toy_network(T=Float64; layout=:volumetric, kwargs...) -> VascularNetwork{T}

Build a small, deterministic [`VascularNetwork`](@ref) shaped like the textbook
microcirculation diagram: a single [`Inlet`](@ref) on the left feeds a bifurcating arterial
tree, its `2^generations` leaves drain through an anastomotic capillary bed, and a venous
tree converges on a single [`Outlet`](@ref) on the right.

Unlike [`generate_network`](@ref) it is pure (no rng) and small, and every count, radius,
and symmetry is known in closed form — a fixture for testing solvers, exporters, and
[`smooth_network`](@ref), which supplies the curvature when composed with it.

Along `x`, a trunk takes the first 15% of `extent`, the arterial fan the next 25%, and the
capillary bed the middle 20%; the venous tree is the mirror image across `x = extent/2`.
Each tree is a binary fan in which every node sits at the centroid of its descendant
leaves. The capillary bed is a ladder — arterial leaf `k` feeds venous leaves `k` and
`k+1` — so it closes loops through the trees. Because every rung leans the same way, the
network is symmetric under a half-turn about the bed center rather than a pure `x`-mirror.

Radii follow Murray's law (`r_parent³ = Σ r_child³`) up from `capillary_radius`: the
vessel entering a generation-`g` node has radius
`capillary_radius * 2^((generations + 1 - g)/3)`. Murray's law holds exactly at every tree
junction; the leaves at the two ends of each ladder feed only one capillary instead of two.

The inlet (velocity `1e-3`) and outlet (pressure `2000`) are degree-1 vertices; every other
vertex is `Interior` with a zero boundary value.

### Layouts
With `n = 2^generations` leaves per tree:

- `:volumetric` (default) — generations alternate between splitting in `y` and in `z`, so
  the leaves form a `2^⌈generations/2⌉ × 2^⌊generations/2⌋` grid in the `y`–`z` plane.
  The leaf order is a Z-order curve, so a few ladder rungs jump across the grid. `4n`
  vertices and `6n - 3` elements, of which `2n - 1` are capillaries.
- `:flat` — the textbook picture: every node in the `z = 0` plane, leaves evenly spaced in
  `y` at pitch `width/n`. Counts match `:volumetric`.
- `:radial` — `copies` of `:flat` rotated about the `x` axis onto planes `π/copies` apart,
  sharing one inlet, one outlet, and both trunks, whose radius grows by `cbrt(copies)` to
  keep Murray's law at the hubs. `2 * (2 + copies * (2n - 2))` vertices and
  `2 * (1 + copies * (2n - 2)) + copies * (2n - 1)` elements.

### Keyword Arguments
- `layout`: one of `:volumetric`, `:flat`, or `:radial` (default `:volumetric`).
- `generations`: bifurcations per tree, `≥ 1`; each tree has `2^generations` leaves
  (default `3`).
- `extent`: inlet-to-outlet distance along `x` (default `200`).
- `width`: transverse spread of the leaves (default `100`).
- `capillary_radius`: radius of every capillary, matching
  `NetworkConfig().radius_default_small_vessels` by default (`2.6`).
- `copies`: number of rotated copies, `≥ 1`; `:radial` only (default `4`).

### Examples
```julia
using Random
net = toy_network()                                   # 32 vertices, 45 elements
net = toy_network(; layout=:flat, generations=2)
net = toy_network(Float32; layout=:radial, copies=3)
curved = smooth_network(MersenneTwister(0), toy_network(; generations=4))
```
"""
function toy_network(::Type{T}=Float64;
    layout=:volumetric,
    generations::Integer=3,
    extent::Real=200,
    width::Real=100,
    capillary_radius::Real=2.6,
    copies::Integer=4,
) where {T<:Real}
    layout ∈ _TOY_LAYOUTS ||
        throw(ArgumentError("layout must be one of $(_TOY_LAYOUTS), got $(repr(layout))"))
    generations ≥ 1 || throw(ArgumentError("generations must be ≥ 1, got $generations"))
    copies ≥ 1 || throw(ArgumentError("copies must be ≥ 1, got $copies"))
    extent > 0 || throw(ArgumentError("extent must be positive, got $extent"))
    width > 0 || throw(ArgumentError("width must be positive, got $width"))
    capillary_radius > 0 ||
        throw(ArgumentError("capillary_radius must be positive, got $capillary_radius"))

    (E, W, rc) = (T(extent), T(width), T(capillary_radius))
    G = generations
    fans = map(axes -> _fan_levels(axes, W), _fan_axes(layout, T, G, copies))
    # Murray's law upward from the capillaries: r³ doubles at each generation toward the hub
    radii = map(g -> rc * T(2)^(T(G + 1 - g) / 3), 0:G)
    trunk_radius = first(radii) * cbrt(T(length(fans)))
    # trunk 0.15, fan 0.25, capillary bed 0.2 of the extent; the venous side mirrors it
    xs = map(g -> E * (T(0.15) + T(0.25) * g / G), 0:G)

    net = VascularNetwork{T}()
    (inlet, arterial) = _add_tree!(net, fans, zero(T), xs, radii, trunk_radius)
    (outlet, venous) = _add_tree!(net, fans, E, map(x -> E - x, xs), radii, trunk_radius)
    set_bc!(net, inlet, Inlet, 1e-3)
    set_bc!(net, outlet, Outlet, 2000)
    foreach((a, v) -> _add_ladder!(net, a, v, rc), arterial, venous)
    return net
end
