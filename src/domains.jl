#--------------------------------------------------------------------------------# Domains

"""
    AbstractDomain

Supertype for the bounding region a network is generated in. The concrete domain a
[`NetworkConfig`](@ref) carries selects which construction [`generate_network`](@ref)
runs.

Subtypes: [`Cube`](@ref), [`HalfSphere`](@ref).
"""
abstract type AbstractDomain end

"""
    Cube()

Cube/box domain: the capillary bed is a Voronoi graph over a square footprint of side
`region_size`, extruded into uniform stacked layers. The default domain.

### Examples
```julia
generate_network(NetworkConfig(domain = Cube()))
```
"""
struct Cube <: AbstractDomain end

"""
    HalfSphere()

Solid-hemisphere domain `{x² + y² + z² ≤ R², z ≥ 0}` of radius `R = sphere_radius`: the
capillary bed fills a stack of shrinking disks from the flat equatorial face (`z = 0`,
radius `R`) up to the apex, and the feeding/draining vessels enter through the curved
dome surface. On this domain `num_layers` acts as an upper bound — the effective layer
count is limited by the geometry (`⌊R / dz⌋ + 1`).

### Examples
```julia
generate_network(NetworkConfig(domain = HalfSphere(), sphere_radius = 120.0))
```
"""
struct HalfSphere <: AbstractDomain end

#--------------------------------------------------------------------------------# Geometry helpers

# Radius of the dome's circular cross-section at height z, i.e. √(R² − z²) clamped to 0
# above the apex. Used to clip each extrusion layer of the half-sphere.
dome_radius(R::T, z::Real) where {T} = sqrt(max(R^2 - T(z)^2, zero(T)))

# Hexagonal-lattice pitch that tiles `area` with about `num_fibers` Voronoi cells. The
# cube derives its pitch from a square (`area = region_size²`); the half-sphere passes
# the disk area `πR²` so `num_fibers` keeps its meaning of "number of capillary cells".
hex_pitch_for_area(area::T, num_fibers::Integer) where {T} =
    sqrt(area) / T(sqrt(3 / 2 * sqrt(3) * num_fibers))

# Orthonormal frame on the sphere at polar angle θ (from +z) and azimuth φ: outward radial
# normal n̂ plus azimuthal (t̂₁) and polar (t̂₂) tangents. Used to place half-sphere I/O
# trunks that pierce the curved surface along n̂ and fan their branches along the tangents.
function dome_frame(::Type{T}, θ::Real, φ::Real) where {T}
    sθ, cθ = sincos(T(θ))
    sφ, cφ = sincos(T(φ))
    n = SVector{3,T}(sθ * cφ, sθ * sφ, cθ)
    t1 = SVector{3,T}(-sφ, cφ, zero(T))
    t2 = SVector{3,T}(cθ * cφ, cθ * sφ, -sθ)
    return n, t1, t2
end
