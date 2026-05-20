[![CI](https://github.com/DerangedIons/IdealizedMicrovasculature.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/DerangedIons/IdealizedMicrovasculature.jl/actions/workflows/CI.yml)
[![Docs Build](https://github.com/DerangedIons/IdealizedMicrovasculature.jl/actions/workflows/Docs.yml/badge.svg)](https://github.com/DerangedIons/IdealizedMicrovasculature.jl/actions/workflows/Docs.yml)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue)](https://DerangedIons.github.io/IdealizedMicrovasculature.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue)](https://DerangedIons.github.io/IdealizedMicrovasculature.jl/dev/)

# IdealizedMicrovasculature.jl

Generate synthetic, idealized myocardial microvascular networks.

The network is built in four stages:

1. **Planar capillary graph** — a jittered hexagonal point set is tessellated and the
   in-region Voronoi edges form a 2D graph.
2. **3D extrusion** — the planar graph is stacked into `num_layers` perturbed layers,
   connected vertically and with a fraction of the inter-layer edges dropped to thin the
   mesh.
3. **Boundaries** — the top and bottom layers are tagged as no-flow.
4. **Feeding/draining vessels** — one draining venule and two feeding arterioles
   (binary branching trees) are appended and connected to the nearest capillaries.

The result is a [`VascularNetwork`](src/types.jl): vertices, line elements, per-element
radii, and per-vertex boundary tags and values.

## Usage

```julia
using Random
using IdealizedMicrovasculature

# reproducible network with the default parameters
net = generate_network(MersenneTwister(0))

# tune the generator
net = generate_network(MersenneTwister(0), NetworkConfig(num_fibers = 60, num_layers = 25))

nvertices(net)   # number of nodes
nelements(net)   # number of vessel segments
```

## Visualization

Plotting is provided by a package extension — load a [Makie](https://docs.makie.org)
backend (`GLMakie` for an interactive window, `CairoMakie` for static output):

```julia
using GLMakie
using IdealizedMicrovasculature

net = generate_network()
plotnetwork(net)   # vessels colored by radius; inlet/outlet nodes highlighted
```

Vessels are drawn as line segments colored by radius. Boundary nodes are highlighted:
crimson for arteriole inlets, blue for the venule outlet, orange for the venule inlet, and
gray for no-flow layers.

This is a Julia port of the idealized network generator from
[`timokoch/dumux-cardiac`](https://github.com/timokoch/dumux-cardiac).
