#--------------------------------------------------------------------------------# generate_network

"""
    generate_network([rng], cfg = NetworkConfig()) -> VascularNetwork

Generate a synthetic idealized myocardial microvascular network.

The construction proceeds in four stages:

 1. A jittered hexagonal point set is tessellated and the in-region Voronoi edges form a
    planar capillary graph.
 2. That graph is extruded into a stacked 3D capillary mesh.
 3. The top and bottom layers are tagged as no-flow boundaries.
 4. One draining venule and two feeding arterioles are appended and connected to the
    capillary mesh; capillaries feeding the trees are widened.

Pass an `rng` for reproducible output.

### Examples
```julia
using Random
net = generate_network(MersenneTwister(0))
net = generate_network(NetworkConfig(num_fibers = 60))
```
"""
function generate_network(rng::AbstractRNG, cfg::NetworkConfig{T}=NetworkConfig()) where {T}
    points, a = hex_grid_points(rng, cfg)
    vertices2d, edges2d = voronoi_network(points, cfg.region_size)

    net = VascularNetwork{T}()
    nvpl = extrude_capillary!(net, rng, vertices2d, edges2d, a, cfg)

    top_offset = (cfg.num_layers - 1) * nvpl
    for i in 1:nvpl
        set_bc!(net, i, NoFlow, 0)
        set_bc!(net, top_offset + i, NoFlow, 0)
    end

    connecting = Int[]
    add_venule!(net, rng, cfg, a, connecting)

    dz = cfg.layer_spacing_fraction * cfg.region_size
    span = cfg.num_layers * dz - T(5)
    sf = cfg.scale_factor_arteriole

    add_arteriole!(net, rng, cfg, connecting;
        start=2 .* SVector{3,T}(70, 115, 0.5 * 0.2 * span),
        finish=2 .* SVector{3,T}(70, 105, 0.5 * 0.3 * span),
        trunk_radius=T(13) * sf,
        o1=2 .* SVector{3,T}(0, -2, 20),
        o2=2 .* SVector{3,T}(-15, -2, 20),
        o11=2 .* SVector{3,T}(0, -2, 20),
        o12=2 .* SVector{3,T}(15, -2, 10),
        o21=2 .* SVector{3,T}(-15, -2, 20),
        o22=2 .* SVector{3,T}(-25, -2, -30),
        exclude_count=8)

    add_arteriole!(net, rng, cfg, connecting;
        start=2 .* SVector{3,T}(-15, 70, 0.5 * 0.8 * span),
        finish=2 .* SVector{3,T}(-5, 70, 0.5 * 0.7 * span),
        trunk_radius=T(12) * sf,
        o1=2 .* SVector{3,T}(-2, -10, -20),
        o2=2 .* SVector{3,T}(-2, 0, -20),
        o11=2 .* SVector{3,T}(-2, -10, -20),
        o12=2 .* SVector{3,T}(2, -30, 30),
        o21=2 .* SVector{3,T}(-2, 0, -20),
        o22=2 .* SVector{3,T}(-2, 20, -10),
        exclude_count=16)

    enforce_min_radius!(net, connecting, cfg)
    return net
end

generate_network(cfg::NetworkConfig=NetworkConfig()) = generate_network(default_rng(), cfg)
