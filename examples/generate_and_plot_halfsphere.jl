#!/usr/bin/env julia
#
# Generate an idealized myocardial microvascular network on the half-sphere
# (dome) domain and visualize it interactively with GLMakie.
#
# Run from the repository root (or anywhere):
#
#     julia examples/generate_and_plot_halfsphere.jl
#
# The script uses its own environment in examples/ and develops the
# IdealizedMicrovasculature package from this repository, so no prior setup is
# needed. The first run precompiles GLMakie and may take a few minutes.

using Pkg
Pkg.activate(@__DIR__)

using Random
using GLMakie
using IdealizedMicrovasculature

# The capillary bed fills a solid hemisphere of radius `sphere_radius`; the
# feeding/draining vessels pierce the curved dome surface. Tune the radius and
# density through the config, e.g. NetworkConfig(domain = HalfSphere(),
# sphere_radius = 120.0, num_fibers = 60).
net = generate_network(MersenneTwister(0), NetworkConfig(domain = HalfSphere()))
@info "Generated network" vertices = nvertices(net) elements = nelements(net)

# Vessels colored by radius; arteriole inlets (crimson), venule outlet (blue),
# venule inlet (orange) and no-flow surfaces (gray) are highlighted.
fig = plotnetwork(net)

# Save a static snapshot next to this script.
snapshot = joinpath(@__DIR__, "idealized_microvasculature_halfsphere.png")
save(snapshot, fig; px_per_unit = 2)
@info "Saved snapshot" path = snapshot

# Open an interactive window — rotate with the mouse, zoom with the scroll
# wheel. Closing the window ends the script.
screen = display(fig)
isinteractive() || wait(screen)
