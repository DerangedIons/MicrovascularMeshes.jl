#!/usr/bin/env julia
#
# Standalone demo of smooth_network: turning the straight-segment output of
# generate_network into organically curved vessels.
#
# Run it directly from a fresh checkout — no environment setup required:
#
#     julia examples/demo_smoothing.jl
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

net = generate_network(MersenneTwister(0))
curved = smooth_network(MersenneTwister(0), net)

println("Straight network: $(nvertices(net)) vertices, $(nelements(net)) elements")
println("Smoothed network: $(nvertices(curved)) vertices, $(nelements(curved)) elements")
println("  growth = $(round(nelements(curved) / nelements(net); digits = 2))x — every element",
    " is subdivided along a cubic curve; junctions, tags, and radii are unchanged")

#--------------------------------------------------------------------------------# Visualize

# Same vertices, only the elements lying entirely inside the box — a close-up crop
# where the curvature is actually visible at print scale.
function crop(network::VascularNetwork{T}, lo, hi) where {T}
    inbox(v) = all(lo .<= v) && all(v .<= hi)
    ks = [k for (k, e) in enumerate(network.elements)
          if inbox(network.vertices[e[1]]) && inbox(network.vertices[e[2]])]
    return VascularNetwork{T}(network.vertices, network.elements[ks], network.radii[ks],
        network.node_tags, network.bcvalues)
end

mins = reduce((a, b) -> min.(a, b), net.vertices)
maxs = reduce((a, b) -> max.(a, b), net.vertices)
center = (mins + maxs) / 2
half = (maxs - mins) * 0.15
sub_lo, sub_hi = center - half, center + half

function panel!(fig, pos, network, title; linewidth = 1.5, show_nodes = true)
    ax = Axis3(fig[pos...]; aspect = :data, title = title,
        azimuth = 0.30π, elevation = π / 8, xlabel = "x", ylabel = "y", zlabel = "z")
    plotnetwork!(ax, network; linewidth = linewidth, show_nodes = show_nodes)
    return ax
end

fig = Figure(size = (1400, 1200))
Label(fig[0, 1:2], "smooth_network — straight segments vs organic curved vessels";
    fontsize = 22, font = :bold)

panel!(fig, (1, 1), net, "Straight (generate_network output)")
panel!(fig, (1, 2), curved, "Smoothed (smooth_network output)")
panel!(fig, (2, 1), crop(net, sub_lo, sub_hi), "Straight — close-up";
    linewidth = 3, show_nodes = false)
panel!(fig, (2, 2), crop(curved, sub_lo, sub_hi), "Smoothed — close-up";
    linewidth = 3, show_nodes = false)

snapshot = joinpath(@__DIR__, "demo_smoothing.png")
save(snapshot, fig; px_per_unit = 2)
println("\nSaved visualization to $snapshot")
