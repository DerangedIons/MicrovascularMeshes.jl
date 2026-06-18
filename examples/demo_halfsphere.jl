#!/usr/bin/env julia
#
# Standalone demo of the half-sphere (dome) domain.
#
# Run it directly from a fresh checkout — no environment setup required:
#
#     julia examples/demo_halfsphere.jl
#
# The script bootstraps its own throwaway environment, develops the
# IdealizedMicrovasculature package from this repository, and adds CairoMakie
# (a headless, no-GPU backend that writes a static PNG). Swap CairoMakie for
# GLMakie below for an interactive, rotatable window. Packages are cached
# globally, so only the first run is slow.

using Pkg

repo_root = dirname(@__DIR__)
Pkg.activate(; temp = true)
Pkg.develop(; path = repo_root)
Pkg.add("CairoMakie")

using Random
using CairoMakie
using IdealizedMicrovasculature

#--------------------------------------------------------------------------------# Generate

R = 100.0
net = generate_network(MersenneTwister(0), NetworkConfig(domain = HalfSphere(), sphere_radius = R))
cube = generate_network(MersenneTwister(0))   # default Cube domain, for comparison

#--------------------------------------------------------------------------------# Summarize

norm3(v) = sqrt(v[1]^2 + v[2]^2 + v[3]^2)
xyradius(v) = sqrt(v[1]^2 + v[2]^2)

println("Half-sphere network (radius R = $R)")
println("  vertices = $(nvertices(net)), elements = $(nelements(net))")

zs = [v[3] for v in net.vertices]
println("  z range  = $(round(minimum(zs); digits = 1)) … $(round(maximum(zs); digits = 1))  (apex at z = R)")
println("  max |p|  = $(round(maximum(norm3, net.vertices); digits = 1))  (I/O trunks start just outside R)")

println("\n  dome cross-section shrinks toward the apex:")
for (lo, hi) in [(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]
    band = [xyradius(v) for v in net.vertices if lo <= v[3] < hi]
    isempty(band) && continue
    println("    z $(lpad(lo, 2))–$(hi): max xy-radius = $(round(maximum(band); digits = 1))")
end

println("\n  feeding/draining vessels pierce the curved surface:")
for (tag, name) in [(Inlet, "arteriole inlet"), (Outlet, "venule outlet"), (InletVenule, "venule inlet")]
    for i in findall(==(tag), net.node_tags)
        v = net.vertices[i]
        println("    $(rpad(name, 16)) |p| = $(round(norm3(v); digits = 1)), z = $(round(v[3]; digits = 1))")
    end
end
println("  no-flow boundary nodes = $(count(==(NoFlow), net.node_tags))  (flat face + dome shell)")

#--------------------------------------------------------------------------------# Visualize

# One panel per camera angle. Vessels are colored by radius; arteriole inlets are
# crimson, the venule outlet blue, the venule inlet orange, no-flow nodes gray.
function panel!(fig, pos, network, title; azimuth, elevation)
    ax = Axis3(fig[pos...]; aspect = :data, title = title,
        azimuth = azimuth, elevation = elevation,
        xlabel = "x", ylabel = "y", zlabel = "z")
    plotnetwork!(ax, network)
    return ax
end

fig = Figure(size = (1400, 1200))
Label(fig[0, 1:2], "Idealized microvasculature — half-sphere domain (R = $R)";
    fontsize = 22, font = :bold)

panel!(fig, (1, 1), net, "Half-sphere — oblique"; azimuth = 0.30π, elevation = π / 8)
panel!(fig, (1, 2), net, "Half-sphere — profile (dome silhouette)"; azimuth = 0.0, elevation = 0.04π)
panel!(fig, (2, 1), net, "Half-sphere — top-down (disk footprint)"; azimuth = 0.30π, elevation = 0.46π)
panel!(fig, (2, 2), cube, "Cube — oblique (for comparison)"; azimuth = 0.30π, elevation = π / 8)

snapshot = joinpath(@__DIR__, "demo_halfsphere.png")
save(snapshot, fig; px_per_unit = 2)
println("\nSaved visualization to $snapshot")
