# IdealizedMicrovasculature.jl

Generates synthetic, idealized myocardial microvascular networks.

## Key files
- module and exports: `src/IdealizedMicrovasculature.jl`
- types: `src/types.jl`
- generation pipeline, in stage order: `src/voronoi.jl` (planar capillary graph) → `src/trees.jl` (arterial/venous trees) → `src/extrude.jl` (2D → 3D) → `src/prune.jl`
- top-level driver: `src/generate.jl`
- region definitions: `src/domains.jl`
- plotting: `ext/IdealizedMicrovasculatureMakieExt.jl` (loads with Makie)

## Gotchas
- The four generation stages are order-dependent and each consumes the previous stage's output — read `src/generate.jl` first to see the pipeline before changing any single stage.
- Plotting lives in a Makie extension, so it is unavailable until the user loads Makie; don't add plotting calls to `src/`.
