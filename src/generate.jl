#--------------------------------------------------------------------------------# generate_network

"""
    generate_network([rng], cfg = NetworkConfig()) -> VascularNetwork

Generate a synthetic idealized myocardial microvascular network.

The bounding region is selected by `cfg.domain` (an [`AbstractDomain`](@ref)). Both domains
share the same construction in four stages:

 1. A jittered hexagonal point set is tessellated and the in-region Voronoi edges form a
    planar capillary graph.
 2. That graph is extruded into a stacked 3D capillary mesh.
 3. The bounding surfaces are tagged as no-flow boundaries.
 4. Feeding arteriole(s) and a draining venule are appended and connected to the capillary
    mesh; capillaries feeding the trees are widened.
 5. Capillaries left disconnected from the main network are pruned (see
    [`prune_disconnected`](@ref)).

On [`Cube`](@ref) the footprint is a square and the mesh is a uniform slab; on
[`HalfSphere`](@ref) the footprint is a disk and the mesh fills a dome whose feeding
vessels pierce the curved surface.

Pass an `rng` for reproducible output.

### Examples
```julia
using Random
net = generate_network(MersenneTwister(0))
net = generate_network(NetworkConfig(num_fibers = 60))
net = generate_network(NetworkConfig(domain = HalfSphere(), sphere_radius = 120.0))
```
"""
generate_network(rng::AbstractRNG, cfg::NetworkConfig=NetworkConfig()) =
    _build_network(cfg.domain, rng, cfg)

generate_network(cfg::NetworkConfig=NetworkConfig()) = generate_network(default_rng(), cfg)

#--------------------------------------------------------------------------------# Cube

function _build_network(::Cube, rng::AbstractRNG, cfg::NetworkConfig{T}) where {T}
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
    return prune_disconnected(net)
end

#--------------------------------------------------------------------------------# HalfSphere

function _build_network(::HalfSphere, rng::AbstractRNG, cfg::NetworkConfig{T}) where {T}
    R = cfg.sphere_radius
    points, a = hex_grid_disk_points(rng, cfg)
    vertices2d, edges2d = voronoi_network(points, R; in_region = c -> c[1]^2 + c[2]^2 <= R^2)

    net = VascularNetwork{T}()
    gid, nlayers = extrude_dome!(net, rng, vertices2d, edges2d, a, R, cfg)

    # No-flow casing: the flat equatorial face (ℓ = 1) plus the curved dome-surface shell,
    # i.e. vertices absent from the next layer up (the rim of their layer) and the apex cap.
    for i in axes(gid, 1), ℓ in 1:nlayers
        g = gid[i, ℓ]
        g == 0 && continue
        (ℓ == 1 || ℓ == nlayers || gid[i, ℓ+1] == 0) && set_bc!(net, g, NoFlow, 0)
    end

    connecting = Int[]
    sf = cfg.scale_factor_arteriole
    δ, ε = T(0.05) * R, T(0.1) * R

    # One feeding arteriole and one draining venule, both piercing the curved surface near
    # the apex: each 15° off the +z axis, on opposite azimuths (30° apart). Trunks start
    # just outside along the outward normal n̂(θ,φ) and finish just inside; the arteriole's
    # bifurcation offsets fan the tree inward (−n̂) and laterally (t̂₁, t̂₂).
    nv, _, _ = dome_frame(T, deg2rad(15), deg2rad(180))
    add_venule!(net, rng, cfg, a, connecting;
        start=(R + δ) * nv, finish=(R - ε) * nv, tag_distal_inlet=false)

    n1, t11, t12 = dome_frame(T, deg2rad(15), deg2rad(0))
    add_arteriole!(net, rng, cfg, connecting;
        start=(R + δ) * n1, finish=(R - ε) * n1,
        trunk_radius=T(13) * sf,
        o1=R * (T(-0.10) * n1 + T(0.10) * t11),
        o2=R * (T(-0.10) * n1 - T(0.10) * t11),
        o11=R * (T(-0.10) * n1 + T(0.08) * t12),
        o12=R * (T(-0.10) * n1 + T(0.12) * t11),
        o21=R * (T(-0.10) * n1 - T(0.08) * t12),
        o22=R * (T(-0.10) * n1 - T(0.12) * t11),
        exclude_count=8)

    enforce_min_radius!(net, connecting, cfg)
    return prune_disconnected(net)
end
