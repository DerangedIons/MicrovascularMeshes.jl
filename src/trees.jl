#--------------------------------------------------------------------------------# Helpers

# Assign a boundary tag and value to a vertex.
function set_bc!(net::VascularNetwork, vertex::Integer, tag::NodeTag, value::Real)
    net.node_tags[vertex] = tag
    net.bcvalues[vertex] = value
    return nothing
end

# Index of the vertex in `net.vertices[1:search_upto]` closest to `point`.
function nearest_vertex(net::VascularNetwork{T}, point::SVector{3,<:Real},
                        search_upto::Integer) where {T}
    best_index = 0
    best_distance = convert(T, Inf)
    @inbounds for i in 1:search_upto
        d = norm(net.vertices[i] - point)
        if d < best_distance
            best_distance = d
            best_index = i
        end
    end
    return best_index
end

#--------------------------------------------------------------------------------# Venule

"""
    add_venule!(net, rng, cfg, a, connecting; start, finish, tag_distal_inlet)

Append the draining venule: a single chain of 11 vertices with a tapering radius, plus
cross-links from every second vertex to the nearest capillary vertex.

The first vertex is tagged [`Outlet`](@ref) (pressure boundary condition). The capillary
vertices the cross-links attach to are appended to `connecting`.

### Keyword Arguments
- `start`, `finish`: trunk endpoints. Default to the `Cube` domain's placement; the
  [`HalfSphere`](@ref) domain passes endpoints that pierce the curved surface.
- `tag_distal_inlet`: when `true` (default), the deep end is tagged [`InletVenule`](@ref)
  (a prescribed venous inflow). Set `false` for a single-outlet venule, leaving the deep
  end an `Interior` capillary-end stub.
"""
function add_venule!(net::VascularNetwork{T}, rng::AbstractRNG, cfg::NetworkConfig{T},
                     a::T, connecting::Vector{Int};
                     start::SVector{3,T}=2 .* SVector{3,T}(0, -5, 0),
                     finish::SVector{3,T}=2 .* SVector{3,T}(100, -5,
                         0.5 * cfg.num_layers * (cfg.layer_spacing_fraction * cfg.region_size) - 5),
                     tag_distal_inlet::Bool=true) where {T}
    nverts = 11
    base = nvertices(net) + 1
    capillary_count = nvertices(net)
    σ = T(0.1) * a
    for i in 0:nverts-1
        p = start + (finish - start) * T(i / nverts)
        p += σ .* SVector{3,T}(_randn(rng, T), _randn(rng, T), _randn(rng, T))
        add_vertex!(net, p)
    end

    set_bc!(net, base, Outlet, 2000)
    tag_distal_inlet && set_bc!(net, base + nverts - 1, InletVenule, 1e-3)

    sf = cfg.scale_factor_venule
    venule_radius = T(20) * sf
    for i in 0:nverts-2
        add_element!(net, base + i, base + i + 1, venule_radius + T(0.1) * _randn(rng, T))
        iseven(i) && (venule_radius -= T(1.5) * sf)
    end

    for i in 1:2:nverts-1
        vertex = base + i
        closest = nearest_vertex(net, net.vertices[vertex], capillary_count)
        push!(connecting, closest)
        add_element!(net, vertex, closest, T(7) * sf + T(0.5) * _randn(rng, T))
    end
    return nothing
end

#--------------------------------------------------------------------------------# Arteriole

"""
    add_arteriole!(net, rng, cfg, connecting; kwargs...)

Append a feeding arteriole: a trunk that bifurcates twice into a binary tree of 8
vertices, whose 4 leaves are each connected to the nearest existing vertex outside the
arteriole trees.

The trunk start is tagged [`Inlet`](@ref) (velocity boundary condition). The vertices the
leaves attach to are appended to `connecting`.

### Keyword Arguments
- `start`, `finish`: trunk endpoints.
- `trunk_radius`: radius of the trunk element.
- `o1`, `o2`: offsets from `finish` to the first-generation daughters.
- `o11`, `o12`: offsets from the first daughter to its two children.
- `o21`, `o22`: offsets from the second daughter to its two children.
- `exclude_count`: number of trailing vertices excluded from the nearest-vertex search
  (the arteriole's own vertices, and those of any earlier arteriole).
"""
function add_arteriole!(net::VascularNetwork{T}, rng::AbstractRNG, cfg::NetworkConfig{T},
                        connecting::Vector{Int};
                        start::SVector{3,T}, finish::SVector{3,T},
                        trunk_radius::Real,
                        o1::SVector{3,T}, o2::SVector{3,T},
                        o11::SVector{3,T}, o12::SVector{3,T},
                        o21::SVector{3,T}, o22::SVector{3,T},
                        exclude_count::Integer) where {T}
    sf = cfg.scale_factor_arteriole
    base = nvertices(net) + 1

    add_vertex!(net, start)
    add_vertex!(net, finish)
    add_element!(net, base, base + 1, trunk_radius)

    daughter1 = finish + o1
    daughter2 = finish + o2
    add_vertex!(net, daughter1)
    add_vertex!(net, daughter2)
    add_element!(net, base + 1, base + 2, T(10) * sf)
    add_element!(net, base + 1, base + 3, T(10) * sf)

    daughter11 = daughter1 + o11
    daughter12 = daughter1 + o12
    add_vertex!(net, daughter11)
    add_vertex!(net, daughter12)
    add_element!(net, base + 2, base + 4, T(7) * sf)
    add_element!(net, base + 2, base + 5, T(7) * sf)

    daughter21 = daughter2 + o21
    daughter22 = daughter2 + o22
    add_vertex!(net, daughter21)
    add_vertex!(net, daughter22)
    add_element!(net, base + 3, base + 6, T(7) * sf)
    add_element!(net, base + 3, base + 7, T(7) * sf)

    leaves = (base + 4, base + 5, base + 6, base + 7)
    probes = (daughter11 + T(0.3) * (daughter11 - daughter1),
              daughter12 + T(0.3) * (daughter12 - daughter1),
              daughter21 + T(0.3) * (daughter21 - daughter2),
              daughter22 + T(0.3) * (daughter22 - daughter2))
    search_upto = nvertices(net) - exclude_count
    for (leaf, probe) in zip(leaves, probes)
        closest = nearest_vertex(net, probe, search_upto)
        push!(connecting, closest)
        add_element!(net, leaf, closest, T(7) * sf)
    end

    set_bc!(net, base, Inlet, 1e-3)
    return nothing
end

#--------------------------------------------------------------------------------# Radius widening

"""
    enforce_min_radius!(net, connecting, cfg)

Widen every element incident to a vertex in `connecting` to at least
`5 * scale_factor_arteriole`, so capillaries feeding the arteriole and venule trees are not
unrealistically thin.
"""
function enforce_min_radius!(net::VascularNetwork{T}, connecting::Vector{Int},
                             cfg::NetworkConfig{T}) where {T}
    rmin = T(5) * cfg.scale_factor_arteriole
    for vertex in connecting
        for (e, element) in enumerate(net.elements)
            if (element[1] == vertex || element[2] == vertex) && net.radii[e] < rmin
                net.radii[e] = rmin
            end
        end
    end
    return nothing
end
