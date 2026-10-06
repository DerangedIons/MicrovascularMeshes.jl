# MicrovascularMeshes.jl

Generates synthetic, idealized myocardial microvascular meshes (geometry, topology, radii, boundary tags) and coordinate-free scaling-law arteriolar trees — geometry only, no flow or transport; the physics lives downstream (Capillarium.jl, MottlingWithaT.jl).

## Key files
- module and exports: `src/MicrovascularMeshes.jl`
- types: `src/types.jl`
- generation pipeline, in stage order: `src/voronoi.jl` (planar capillary graph) → `src/trees.jl` (arterial/venous trees) → `src/extrude.jl` (2D → 3D) → `src/prune.jl`
- top-level driver: `src/generate.jl`
- hand-built deterministic test/demo fixture: `src/toy.jl` (`toy_network`, independent of the generation pipeline)
- region definitions: `src/domains.jl`
- scaling-law arteriolar trees: `src/scaling.jl` (`scaling_tree`, coordinate-free, independent of the mesh pipeline; its docstring carries the literature)
- plotting: `ext/MicrovascularMeshesMakieExt.jl` (loads with Makie)

## Gotchas
- The four generation stages are order-dependent and each consumes the previous stage's output — read `src/generate.jl` first to see the pipeline before changing any single stage.
- Plotting lives in a Makie extension, so it is unavailable until the user loads Makie; don't add plotting calls to `src/`.
