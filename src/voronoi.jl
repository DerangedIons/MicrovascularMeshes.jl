#--------------------------------------------------------------------------------# Hex grid

# Half-open float range matching numpy's `arange(start, stop, step)`.
function _arange(start::Real, stop::Real, step::Real)
    n = max(0, ceil(Int, (stop - start) / step))
    return [start + k * step for k in 0:n-1]
end

# Standard-normal draw converted to the element type `T`. One idiom for every `T`:
# `randn(rng, T)` is undefined for non-`AbstractFloat` reals, `convert` is not.
_randn(rng::AbstractRNG, ::Type{T}) where {T<:Real} = convert(T, randn(rng))

# Uniform-[0,1) draw converted to the element type `T`, for the same reason.
_rand(rng::AbstractRNG, ::Type{T}) where {T<:Real} = convert(T, rand(rng))

# Lay two interleaved rectangular grids (a hexagonal lattice) over the coordinate lines
# `x`,`y` with pitch `a`, then jitter every point by Gaussian noise of stddev `σ`. The
# domain-specific seeders below differ only in the coordinate lines and pitch they pass.
function _hex_lattice(rng::AbstractRNG, ::Type{T}, x::Vector{T}, y::Vector{T},
                      a::T, σ::T) where {T}
    points = NTuple{2,T}[]
    for yi in y, xi in x
        push!(points, (xi, yi))
    end
    for yi in y .+ T(sqrt(3)) * a / 2, xi in x .+ 3a / 2
        push!(points, (xi, yi))
    end
    for k in eachindex(points)
        px, py = points[k]
        points[k] = (px + σ * _randn(rng, T), py + σ * _randn(rng, T))
    end
    return points
end

"""
    hex_grid_points(rng, cfg) -> (points, a)

Build the jittered hexagonal seed points for the Voronoi tessellation over the [`Cube`](@ref)
domain's square footprint.

Two interleaved rectangular grids form a hexagonal lattice covering
`[-a, region_size + a]²` (slightly larger than the region so the tessellation's interior
cells are well defined). Every point is then perturbed by Gaussian noise with standard
deviation `random_perturbation * a`.

Returns the seed `points` as a `Vector{NTuple{2,T}}` and the grid pitch `a`.
"""
function hex_grid_points(rng::AbstractRNG, cfg::NetworkConfig{T}) where {T}
    rs = cfg.region_size
    a = rs / T(sqrt(3 / 2 * sqrt(3) * cfg.num_fibers))
    x = _arange(-a, rs + a, 3a)
    y = _arange(-a, rs + a, T(sqrt(3)) * a)
    σ = cfg.random_perturbation * a
    return _hex_lattice(rng, T, x, y, a, σ), a
end

"""
    hex_grid_disk_points(rng, cfg) -> (points, a)

Build the jittered hexagonal seed points for the [`HalfSphere`](@ref) domain: an
origin-centered hexagonal lattice covering `[-R - a, R + a]²` (`R = sphere_radius`). The
pitch `a` is derived from the disk area `πR²` so `num_fibers` keeps its meaning of "number
of capillary cells". Membership in the disk itself is enforced downstream by
[`voronoi_network`](@ref)'s `in_region` predicate, mirroring the square test on the
`Cube` domain.

Returns the seed `points` as a `Vector{NTuple{2,T}}` and the grid pitch `a`.
"""
function hex_grid_disk_points(rng::AbstractRNG, cfg::NetworkConfig{T}) where {T}
    R = cfg.sphere_radius
    a = hex_pitch_for_area(T(π) * R^2, cfg.num_fibers)
    x = _arange(-R - a, R + a, 3a)
    y = _arange(-R - a, R + a, T(sqrt(3)) * a)
    σ = cfg.random_perturbation * a
    return _hex_lattice(rng, T, x, y, a, σ), a
end

#--------------------------------------------------------------------------------# Voronoi network

"""
    voronoi_network(points, region_size) -> (vertices, edges)

Compute the Voronoi tessellation of `points` and extract the planar graph of finite
Voronoi edges that lie fully inside `[0, region_size]²`.

`DelaunayTriangulation.jl` returns each Voronoi cell as a closed ring of vertex indices,
with negative indices marking unbounded (ghost) cells. Edges touching a ghost vertex are
discarded — the analogue of filtering scipy's `-1` ridge vertices. The surviving edges are
deduplicated (each interior edge is shared by two cells) and their endpoints remapped to
compact 1-based indices.

`triangulate` inserts points in a randomized order, so the tessellation's vertex labels
vary between runs even for identical input. The output is therefore canonicalized: vertices
are sorted by coordinate and edges by endpoint, making the function a pure, reproducible
mapping of `points` to a graph.

`DelaunayTriangulation.jl`'s robust geometric predicates only operate in `Float64`, so
the tessellation runs there regardless of `T`; the extracted vertices are then converted
back to `T`. The discrete graph (edges and vertex ordering) is `T`-independent anyway.

`in_region(c)` decides whether a Voronoi vertex at coordinate `c` is kept; an edge
survives only when both endpoints satisfy it. The default is the `Cube` domain's square
test `[0, region_size]²`; the [`HalfSphere`](@ref) domain passes a disk test instead.

Returns `vertices::Vector{SVector{2,T}}` and `edges::Vector{SVector{2,Int}}`.
"""
function voronoi_network(points::Vector{NTuple{2,T}}, region_size::Real;
                         in_region = c -> 0.0 <= c[1] <= region_size &&
                                          0.0 <= c[2] <= region_size) where {T<:Real}
    tess = DT.voronoi(DT.triangulate(NTuple{2,Float64}.(points)))
    polygon_points = DT.get_polygon_points(tess)

    edge_set = Set{Tuple{Int,Int}}()
    for ring in values(DT.get_polygons(tess))
        for k in 1:length(ring)-1
            p, q = ring[k], ring[k+1]
            (p < 1 || q < 1) && continue
            push!(edge_set, minmax(p, q))
        end
    end

    remap = Dict{Int,Int}()
    vertices = SVector{2,T}[]
    edges = SVector{2,Int}[]
    add_vertex(i) = get!(remap, i) do
        c = polygon_points[i]
        push!(vertices, SVector{2,T}(c[1], c[2]))
        length(vertices)
    end
    for (p, q) in edge_set
        (in_region(polygon_points[p]) && in_region(polygon_points[q])) || continue
        push!(edges, SVector{2,Int}(add_vertex(p), add_vertex(q)))
    end

    # Canonicalize so the result is independent of the triangulation's vertex labels.
    order = sortperm(vertices; by = v -> (v[1], v[2]))
    rank = invperm(order)
    canonical_vertices = vertices[order]
    canonical_edges = [SVector{2,Int}(minmax(rank[e[1]], rank[e[2]])...) for e in edges]
    sort!(canonical_edges; by = e -> (e[1], e[2]))
    return canonical_vertices, canonical_edges
end
