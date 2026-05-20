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
function extrude_capillary!(net::VascularNetwork, rng::AbstractRNG,
                            vertices2d::Vector{SVector{2,Float64}},
                            edges2d::Vector{SVector{2,Int}},
                            a::Real, cfg::NetworkConfig)
    nvpl = length(vertices2d)
    nlayers = cfg.num_layers
    dz = cfg.layer_spacing_fraction * cfg.region_size

    base = [SVector{3,Float64}(v[1], v[2], 0.0) for v in vertices2d]
    for p in base
        add_vertex!(net, p)
    end

    σxy = 0.1 * a
    σz = 0.005 * cfg.region_size
    for layer in 1:nlayers-1
        z_offset = layer * dz
        perturb = layer != nlayers - 1
        for p in base
            q = SVector{3,Float64}(p[1], p[2], p[3] + z_offset)
            if perturb
                q += SVector{3,Float64}(σxy * randn(rng), σxy * randn(rng), σz * randn(rng))
            end
            add_vertex!(net, q)
        end
    end

    μ = log(cfg.radius_default_small_vessels)
    σ = cfg.radius_stddev_small_vessels
    sample_radius() = exp(μ + σ * randn(rng))

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
