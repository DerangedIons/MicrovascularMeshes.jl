#!/usr/bin/env julia
#
# Generate an idealized myocardial microvascular network and visualize it
# interactively with GLMakie.
#
# Run from the repository root (or anywhere):
#
#     julia examples/generate_and_plot.jl
#
# The script uses its own environment in examples/ and develops the
# IdealizedMicrovasculature package from this repository, so no prior setup is
# needed. The first run precompiles GLMakie and may take a few minutes.

using Pkg
Pkg.activate(@__DIR__)

using Random
using GLMakie
using IdealizedMicrovasculature

# Reproducible network. Drop the RNG argument for a different one each run, or
# pass a NetworkConfig to change the size/density, e.g.
#     generate_network(MersenneTwister(0), NetworkConfig(num_fibers = 60))
net = generate_network(MersenneTwister(0))
@info "Generated network" vertices = nvertices(net) elements = nelements(net)

# Replace the straight segments with organically curved vessels.
net = smooth_network(MersenneTwister(0), net)
@info "Smoothed network" vertices = nvertices(net) elements = nelements(net)

# Vessels colored by radius; arteriole inlets (crimson), venule outlet (blue),
# venule inlet (orange) and no-flow layers (gray) are highlighted.
fig = plotnetwork(net)

# Save a static snapshot next to this script.
snapshot = joinpath(@__DIR__, "idealized_microvasculature.png")
save(snapshot, fig; px_per_unit = 2)
@info "Saved snapshot" path = snapshot

# Open an interactive window — rotate with the mouse, zoom with the scroll
# wheel. Closing the window ends the script.
screen = display(fig)
isinteractive() || wait(screen)
