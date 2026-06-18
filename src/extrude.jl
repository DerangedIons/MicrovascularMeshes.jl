#--------------------------------------------------------------------------------# Edge dropping

"""
    drop_horizontal_edges(rng, nedges, drop_fraction) -> kept

Randomly select which of `nedges` inter-layer edges to keep on one layer. A fraction
`drop_fraction` (truncated to an integer count) is dropped; the sorted indices of the
surviving edges are returned.
"""
function drop_horizontal_edges(rng::AbstractRNG, nedges::Integer, drop_fraction::Real)
    ndrop = floor(Int, drop_fraction * nedges)
    return sort!(randperm(rng, nedges)[ndrop+1:end])
end

#--------------------------------------------------------------------------------# Extrusion

"""
    extrude_capillary!(net, rng, vertices2d, edges2d, a, cfg) -> nvertices_per_layer

Extrude the planar Voronoi graph into a 3D capillary mesh, writing vertices and elements
into `net`.

The planar graph is copied into `num_layers` stacked layers spaced
`layer_spacing_fraction * region_size` apart. Every layer except the bottom (`z = 0`) and
the top is jittered. Vertical elements connect each vertex to its copy in the next layer;
"horizontal" elements re-use the planar edges but span consecutive layers, with a fraction
`drop_fraction_hex_grid` dropped per layer to thin the network. Every element is given a
lognormally distributed radius.

Returns the number of vertices per layer.
"""
function extrude_capillary!(net::VascularNetwork{T}, rng::AbstractRNG,
                            vertices2d::Vector{SVector{2,T}},
                            edges2d::Vector{SVector{2,Int}},
                            a::T, cfg::NetworkConfig{T}) where {T}
    nvpl = length(vertices2d)
    nlayers = cfg.num_layers
    dz = cfg.layer_spacing_fraction * cfg.region_size

    base = [SVector{3,T}(v[1], v[2], 0) for v in vertices2d]
    for p in base
        add_vertex!(net, p)
    end

    σxy = T(0.1) * a
    σz = T(0.005) * cfg.region_size
    for layer in 1:nlayers-1
        z_offset = layer * dz
        perturb = layer != nlayers - 1
        for p in base
            q = SVector{3,T}(p[1], p[2], p[3] + z_offset)
            if perturb
                q += SVector{3,T}(σxy * _randn(rng, T), σxy * _randn(rng, T), σz * _randn(rng, T))
            end
            add_vertex!(net, q)
        end
    end

    μ = log(cfg.radius_default_small_vessels)
    σ = cfg.radius_stddev_small_vessels
    sample_radius() = exp(μ + σ * _randn(rng, T))

    # Vertical elements connect a vertex to its copy in the next layer.
    for layer in 1:nlayers-1, i in 1:nvpl
        add_element!(net, i + (layer - 1) * nvpl, i + layer * nvpl, sample_radius())
    end

    # Horizontal elements re-use each planar edge but span consecutive layers.
    for layer in 0:nlayers-2
        for idx in drop_horizontal_edges(rng, length(edges2d), cfg.drop_fraction_hex_grid)
            e = edges2d[idx]
            add_element!(net, e[1] + layer * nvpl, e[2] + (layer + 1) * nvpl, sample_radius())
        end
    end
    return nvpl
end

#--------------------------------------------------------------------------------# Dome extrusion

"""
    extrude_dome!(net, rng, vertices2d, edges2d, a, R, cfg) -> (gid, nlayers)

Extrude the planar disk Voronoi graph into a solid-hemisphere capillary mesh of radius
`R`, writing vertices and elements into `net`.

Layers are stacked at `z = (ℓ-1)·dz` with `dz = layer_spacing_fraction · R`, and each
layer is clipped to the dome's cross-section: planar vertex `i` is kept in layer `ℓ` only
when `‖vertices2d[i]‖ ≤ √(R² − z²)`. Because that radius shrinks with height the active
vertex sets are nested, so the layers form a dome rather than a prism. The layer count is
geometry-limited to `min(num_layers, ⌊R/dz⌋ + 1)` so no layer climbs past the apex.

Connectivity mirrors [`extrude_capillary!`](@ref): vertical elements join a vertex to its
copy in the next layer, horizontal elements re-use each planar edge across consecutive
layers (a fraction dropped per layer) — but an element is added only when both of its
endpoints survive the clip in their respective layers, which also guarantees no vertical
dangles off an absent vertex. Every element is given a lognormally distributed radius.

Returns `gid`, an `nvertices_per_disk × nlayers` matrix mapping `(planar vertex, layer)`
to the global network vertex index (`0` where the vertex is clipped away), and the
effective `nlayers`. `gid` lets the caller identify the flat-face and dome-surface
boundary vertices.
"""
function extrude_dome!(net::VascularNetwork{T}, rng::AbstractRNG,
                       vertices2d::Vector{SVector{2,T}},
                       edges2d::Vector{SVector{2,Int}},
                       a::T, R::T, cfg::NetworkConfig{T}) where {T}
    np = length(vertices2d)
    dz = cfg.layer_spacing_fraction * R
    nlayers = min(cfg.num_layers, floor(Int, R / dz) + 1)

    gid = zeros(Int, np, nlayers)
    σxy = T(0.1) * a
    σz = T(0.005) * R
    for ℓ in 1:nlayers
        z = (ℓ - 1) * dz
        rℓ = dome_radius(R, z)
        perturb = ℓ != 1 && ℓ != nlayers
        for i in 1:np
            v = vertices2d[i]
            (v[1]^2 + v[2]^2 <= rℓ^2) || continue
            p = SVector{3,T}(v[1], v[2], z)
            if perturb
                p += SVector{3,T}(σxy * _randn(rng, T), σxy * _randn(rng, T), σz * _randn(rng, T))
            end
            gid[i, ℓ] = add_vertex!(net, p)
        end
    end

    μ = log(cfg.radius_default_small_vessels)
    σ = cfg.radius_stddev_small_vessels
    sample_radius() = exp(μ + σ * _randn(rng, T))

    # Vertical elements connect a vertex to its copy in the next layer, when present in both.
    for ℓ in 1:nlayers-1, i in 1:np
        (gid[i, ℓ] != 0 && gid[i, ℓ+1] != 0) || continue
        add_element!(net, gid[i, ℓ], gid[i, ℓ+1], sample_radius())
    end

    # Horizontal elements re-use each planar edge across consecutive layers, skipping any
    # whose endpoints were clipped. The per-layer drop draw matches extrude_capillary!.
    for ℓ in 1:nlayers-1
        for idx in drop_horizontal_edges(rng, length(edges2d), cfg.drop_fraction_hex_grid)
            e = edges2d[idx]
            g1, g2 = gid[e[1], ℓ], gid[e[2], ℓ+1]
            (g1 != 0 && g2 != 0) || continue
            add_element!(net, g1, g2, sample_radius())
        end
    end
    return gid, nlayers
end
