#--------------------------------------------------------------------------------# Smoothing

# Unit vector along `v`, or the zero vector when `v` has zero norm (degenerate guard).
function _unit(v::SVector{3,T}) where {T}
    m = norm(v)
    return iszero(m) ? v : v / m
end

# Median of a non-empty vector, without a Statistics dependency.
function _median(xs::Vector{T}) where {T}
    s = sort(xs)
    n = length(s)
    mid = (n + 1) >> 1
    return isodd(n) ? s[mid] : (s[mid] + s[mid+1]) / 2
end

# Rotate the unit vector `t` by a random angle in [0, θmax) about a random axis
# perpendicular to `t`. Exactly two `_rand` draws, so output is rng-deterministic.
function _tilt(rng::AbstractRNG, t::SVector{3,T}, θmax::T) where {T}
    (iszero(θmax) || iszero(t)) && return t
    ax, ay, az = abs(t[1]), abs(t[2]), abs(t[3])
    e = ax <= ay && ax <= az ? SVector{3,T}(1, 0, 0) :
        (ay <= az ? SVector{3,T}(0, 1, 0) : SVector{3,T}(0, 0, 1))
    u = _unit(cross(t, e))
    w = cross(t, u)
    θ = θmax * _rand(rng, T)
    φ = 2 * T(π) * _rand(rng, T)
    return cos(θ) * t + sin(θ) * (cos(φ) * u + sin(φ) * w)
end

# One unit tangent per element end, oriented OUTGOING from that end's vertex:
# `tstart[k]` leaves `elements[k][1]`, `tend[k]` leaves `elements[k][2]`. Assigned
# per vertex so the ends meeting at a vertex can share a tangent direction:
# - degree 1: the element's own chord, randomly tilted.
# - degree 2: one shared tangent through the vertex (never tilted).
# - degree ≥ 3: the most-collinear pair of chords shares a tangent (the "through
#   vessel", never tilted); every other end gets its own chord, randomly tilted.
function _end_tangents(rng::AbstractRNG, net::VascularNetwork{T}, θmax::T) where {T}
    nv, ne = nvertices(net), nelements(net)
    incident = [Int[] for _ in 1:nv]
    for (k, e) in enumerate(net.elements)
        push!(incident[e[1]], k)
        push!(incident[e[2]], k)
    end

    tstart = fill(zero(SVector{3,T}), ne)
    tend = fill(zero(SVector{3,T}), ne)
    setat!(k, v, τ) = net.elements[k][1] == v ? (tstart[k] = τ) : (tend[k] = τ)
    other(k, v) = (e = net.elements[k]; e[1] == v ? e[2] : e[1])
    chord(k, v) = _unit(net.vertices[other(k, v)] - net.vertices[v])

    for v in 1:nv
        elems = incident[v]
        deg = length(elems)
        if deg == 1
            k = elems[1]
            setat!(k, v, _tilt(rng, chord(k, v), θmax))
        elseif deg == 2
            k1, k2 = elems[1], elems[2]
            t = _unit(net.vertices[other(k2, v)] - net.vertices[other(k1, v)])
            setat!(k1, v, -t)
            setat!(k2, v, t)
        elseif deg >= 3
            ks = Int[]
            ds = SVector{3,T}[]
            for k in elems
                d = chord(k, v)
                iszero(d) && continue
                push!(ks, k)
                push!(ds, d)
            end
            # Through pair: outgoing chords closest to opposite (dot nearest -1).
            besti, bestj, bestdot = 0, 0, T(2)
            for p in 1:length(ks)-1, q in p+1:length(ks)
                dpq = dot(ds[p], ds[q])
                dpq < bestdot && ((besti, bestj, bestdot) = (p, q, dpq))
            end
            if besti != 0
                t = _unit(ds[besti] - ds[bestj])
                if iszero(t)  # all incident chords parallel — keep own directions
                    setat!(ks[besti], v, ds[besti])
                    setat!(ks[bestj], v, ds[bestj])
                else
                    setat!(ks[besti], v, t)
                    setat!(ks[bestj], v, -t)
                end
            end
            for p in eachindex(ks)
                (p == besti || p == bestj) && continue
                setat!(ks[p], v, _tilt(rng, ds[p], θmax))
            end
        end
    end
    return tstart, tend
end

"""
    smooth_network(rng, net; max_tilt = π/9, spacing = nothing) -> VascularNetwork
    smooth_network(net; kwargs...) -> VascularNetwork

Return a refined copy of `net` in which every straight element is replaced by a sampled
cubic Hermite curve, giving the vessels an organic, curved shape. Curves are
tangent-continuous through pass-through (degree-2) vertices and along the most-collinear
pair of elements at each junction (the "through vessel"); the remaining branches take off
along their own chord directions. Organic variation comes from small random tilts of the
free end tangents, drawn from `rng`.

The transformation is pure geometric refinement:
- Original vertices keep their indices (`1:nvertices(net)`), positions, tags, and
  boundary values; new sample vertices are tagged `Interior` with zero BC value.
- Every sub-element inherits its parent element's radius.
- Topology is unchanged — no vertex moves and no connection is added or removed, elements
  are only subdivided. Elements shorter than `spacing` stay single straight segments.

### Keyword Arguments
- `max_tilt`: maximum random tilt angle in radians applied to free end tangents; `0`
  disables the organic wiggle (default `π/9` ≈ 20°).
- `spacing`: target sub-element length in world units; each element is split into
  `ceil(length/spacing)` pieces. The default `nothing` uses a quarter of the median
  element length, growing the element count roughly 4-6×.

### Examples
```julia
using Random
net = generate_network(MersenneTwister(0))
curved = smooth_network(MersenneTwister(0), net)
```
"""
function smooth_network(rng::AbstractRNG, net::VascularNetwork{T};
                        max_tilt::Real = π / 9,
                        spacing::Union{Nothing,Real} = nothing) where {T}
    max_tilt >= 0 || throw(ArgumentError("max_tilt must be non-negative, got $max_tilt"))
    spacing === nothing || spacing > 0 ||
        throw(ArgumentError("spacing must be positive, got $spacing"))

    out = VascularNetwork{T}(copy(net.vertices), SVector{2,Int}[], T[],
                             copy(net.node_tags), copy(net.bcvalues))
    nelements(net) == 0 && return out

    lengths = [norm(net.vertices[e[2]] - net.vertices[e[1]]) for e in net.elements]
    h = spacing === nothing ? _median(lengths) / 4 : T(spacing)
    if iszero(h)  # every chord degenerate — nothing to curve
        append!(out.elements, net.elements)
        append!(out.radii, net.radii)
        return out
    end

    tstart, tend = _end_tangents(rng, net, T(max_tilt))

    for (k, e) in enumerate(net.elements)
        i, j = e[1], e[2]
        L = lengths[k]
        n = max(1, ceil(Int, L / h))
        if n == 1  # too short to subdivide (covers L == 0)
            add_element!(out, i, j, net.radii[k])
            continue
        end
        P0, P1 = net.vertices[i], net.vertices[j]
        m0 = L * tstart[k]
        m1 = -L * tend[k]  # tend points backward along the i -> j traversal
        prev = i
        for s in 1:n-1
            t = T(s / n)
            t2, t3 = t * t, t * t * t
            p = (2t3 - 3t2 + 1) * P0 + (t3 - 2t2 + t) * m0 +
                (-2t3 + 3t2) * P1 + (t3 - t2) * m1
            idx = add_vertex!(out, p)
            add_element!(out, prev, idx, net.radii[k])
            prev = idx
        end
        add_element!(out, prev, j, net.radii[k])
    end
    return out
end

smooth_network(net::VascularNetwork; kwargs...) = smooth_network(default_rng(), net; kwargs...)
