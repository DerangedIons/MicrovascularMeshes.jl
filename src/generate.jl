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
function generate_network(rng::AbstractRNG, cfg::NetworkConfig=NetworkConfig())
    points, a = hex_grid_points(rng, cfg)
    vertices2d, edges2d = voronoi_network(points, cfg.region_size)

    net = VascularNetwork()
    nvpl = extrude_capillary!(net, rng, vertices2d, edges2d, a, cfg)

    top_offset = (cfg.num_layers - 1) * nvpl
    for i in 1:nvpl
        set_bc!(net, i, NoFlow, 0.0)
        set_bc!(net, top_offset + i, NoFlow, 0.0)
    end

    connecting = Int[]
    add_venule!(net, rng, cfg, a, connecting)

    dz = cfg.layer_spacing_fraction * cfg.region_size
    span = cfg.num_layers * dz - 5.0
    sf = cfg.scale_factor_arteriole

    add_arteriole!(net, rng, cfg, connecting;
        start=2 .* SVector{3,Float64}(70.0, 115.0, 0.5 * 0.2 * span),
        finish=2 .* SVector{3,Float64}(70.0, 105.0, 0.5 * 0.3 * span),
        trunk_radius=13.0 * sf,
        o1=2 .* SVector{3,Float64}(0.0, -2.0, 20.0),
        o2=2 .* SVector{3,Float64}(-15.0, -2.0, 20.0),
        o11=2 .* SVector{3,Float64}(0.0, -2.0, 20.0),
        o12=2 .* SVector{3,Float64}(15.0, -2.0, 10.0),
        o21=2 .* SVector{3,Float64}(-15.0, -2.0, 20.0),
        o22=2 .* SVector{3,Float64}(-25.0, -2.0, -30.0),
        exclude_count=8)

    add_arteriole!(net, rng, cfg, connecting;
        start=2 .* SVector{3,Float64}(-15.0, 70.0, 0.5 * 0.8 * span),
        finish=2 .* SVector{3,Float64}(-5.0, 70.0, 0.5 * 0.7 * span),
        trunk_radius=12.0 * sf,
        o1=2 .* SVector{3,Float64}(-2.0, -10.0, -20.0),
        o2=2 .* SVector{3,Float64}(-2.0, 0.0, -20.0),
        o11=2 .* SVector{3,Float64}(-2.0, -10.0, -20.0),
        o12=2 .* SVector{3,Float64}(2.0, -30.0, 30.0),
        o21=2 .* SVector{3,Float64}(-2.0, 0.0, -20.0),
        o22=2 .* SVector{3,Float64}(-2.0, 20.0, -10.0),
        exclude_count=16)

    enforce_min_radius!(net, connecting, cfg)
    return net
end

generate_network(cfg::NetworkConfig=NetworkConfig()) = generate_network(default_rng(), cfg)
