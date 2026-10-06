[![CI](https://github.com/DerangedIons/MicrovascularMeshes.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/DerangedIons/MicrovascularMeshes.jl/actions/workflows/CI.yml)
[![Docs Build](https://github.com/DerangedIons/MicrovascularMeshes.jl/actions/workflows/Docs.yml/badge.svg)](https://github.com/DerangedIons/MicrovascularMeshes.jl/actions/workflows/Docs.yml)
[![Stable Docs](https://img.shields.io/badge/docs-stable-blue)](https://DerangedIons.github.io/MicrovascularMeshes.jl/stable/)
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue)](https://DerangedIons.github.io/MicrovascularMeshes.jl/dev/)

# MicrovascularMeshes.jl

Generate synthetic, idealized myocardial microvascular meshes: the vessel geometry, topology, radii, and boundary tags only. No flow or transport is solved here; see [Capillarium.jl](https://github.com/DerangedIons/Capillarium.jl) for that.

It also generates coordinate-free arteriolar trees whose diameters and lengths follow measured scaling laws (`scaling_tree`, see [below](#scaling-law-trees)).

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

The bounding region is set by `cfg.domain`: the default [`Cube`](src/domains.jl) gives the
square-footprint slab above, while [`HalfSphere`](src/domains.jl) fills a solid hemisphere
of radius `sphere_radius` whose feeding/draining vessels pierce the curved dome surface.

## Usage

```julia
using Random
using MicrovascularMeshes

# reproducible network with the default parameters
net = generate_network(MersenneTwister(0))

# tune the generator
net = generate_network(MersenneTwister(0), NetworkConfig(num_fibers = 60, num_layers = 25))

# generate the bed inside a half-sphere instead of a cube
net = generate_network(MersenneTwister(0), NetworkConfig(domain = HalfSphere(), sphere_radius = 120.0))

nvertices(net)   # number of nodes
nelements(net)   # number of vessel segments
```

## Scaling-law trees

`scaling_tree` builds a symmetric bifurcating tree whose diameters follow a junction exponent γ (Murray's 3 by default), whose lengths follow `L ∝ D^β`, and whose vessels each scatter about their generation's mean by a lognormal coefficient of variation. It returns a `ScalingTree` (parents, generation, diameter, length) with no coordinates, for callers that need only the hydraulics. The docstring says which literature each exponent and spread comes from.

```julia
tree = scaling_tree(MersenneTwister(0), 8; leaf_diameter=30.0, leaf_length=600.0, β=1.2,
                    cv_diameter=0.1, cv_length=0.3)
```

## Visualization

Plotting is provided by a package extension — load a [Makie](https://docs.makie.org)
backend (`GLMakie` for an interactive window, `CairoMakie` for static output):

```julia
using GLMakie
using MicrovascularMeshes

net = generate_network()
plotnetwork(net)   # vessels colored by radius; inlet/outlet nodes highlighted
```

Vessels are drawn as line segments colored by radius. Boundary nodes are highlighted:
crimson for arteriole inlets, blue for the venule outlet, orange for the venule inlet, and
gray for no-flow layers.

For a complete runnable example, see [`examples/generate_and_plot.jl`](examples/generate_and_plot.jl) —
it sets up its own environment and opens an interactive 3D window:

```
julia examples/generate_and_plot.jl
```

This is a Julia port of the idealized network generator from
[`timokoch/dumux-cardiac`](https://github.com/timokoch/dumux-cardiac).
