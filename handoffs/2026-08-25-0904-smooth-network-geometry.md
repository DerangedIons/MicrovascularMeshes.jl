---
slug: smooth-network-geometry
created: 2026-08-25-0904
status: done  # implemented as smooth_network in commit 4e2ee48
---

# Handoff: Smoothing function to make network tubes organic/curved

## Goal / why this matters
The generation pipeline produces vessel geometry as straight line segments between nodes, which looks artificial. Kyle wants a new function that takes the output geometry and "smooths" it so tubes are organic-shaped/curved.

## Background & current state
Planning had just started (brainstorming classification: **bounded** — one new exported post-processing function following the existing pipeline-stage pattern). No code was written and **no design was approved** — the brainstorming approval gate was never passed; present a short design and get Kyle's yes before implementing.

Kyle already answered the key requirement question: the smoothed output is for **both** visualization and downstream geometry/mesh export. So smoothing must be a real network transformation — same network type in, same type out, with edges subdivided into polyline sub-elements and radii carried along — NOT a plotting-only trick in the Makie extension.

A structural exploration of the codebase was completed; its findings are baked into the sections below (verified 2026-08-25).

## Key files / locations
- Types: `src/types.jl` — output type is `VascularNetwork{T<:Real}` (`src/types.jl:92-102`) with fields `vertices::Vector{SVector{3,T}}`, `elements::Vector{SVector{2,Int}}` (1-based vertex index pairs), `radii::Vector{T}` (**per element**), `node_tags::Vector{NodeTag}` and `bcvalues::Vector{T}` (**per vertex**). Invariants: `length(radii) == length(elements)`, `length(node_tags) == length(bcvalues) == length(vertices)`.
- Mutation helpers to reuse: `add_vertex!` (`src/types.jl:116`, defaults tag `Interior`, bc `zero(T)`) and `add_element!` (`src/types.jl:129`, pushes element + radius together); `nvertices`/`nelements` (`src/types.jl:105-108`).
- Driver: `src/generate.jl` — entry point is `generate_network` (`src/generate.jl:34-37`); both `_build_network` methods (Cube `:41-87`, HalfSphere `:91-133`) end with `return prune_disconnected(net)`. **The smoothing stage slots at those two call sites** (or stays user-invoked — Kyle to decide).
- Closest pattern: `src/prune.jl:22-55` — `prune_disconnected(net::VascularNetwork{T}) -> VascularNetwork{T}`, pure/out-of-place rebuild with an old→new index remap. Note its early-outs return the *same object* (aliasing, documented at `src/prune.jl:20`) — decide explicitly whether smooth mirrors that.
- Node tags: `@enum NodeTag` = Interior/NoFlow/Inlet/Outlet/InletVenule (`src/types.jl:15-21`). **There is no arterial/venous/capillary label and no per-element vessel-type label** — radius is the only post-hoc signal to treat big vessels differently from capillaries.
- Plotting: `ext/MicrovascularMeshesMakieExt.jl` — `plotnetwork!` draws flat 2-point `linesegments!` with the element radius as a *color* channel (`:21-30`); no tubes/meshes. A subdivided network renders unchanged (just more, shorter segments). Function stubs live in `src/MicrovascularMeshes.jl:43-44`.
- Tests: `test/runtests.jl` includes one file per module inside one outer `@testset`, with `const MM = MicrovascularMeshes`, `TEST_TYPES = (Float64, Float32)`, and an `ncomponents` helper (`test/runtests.jl:9-37`). `test/test_prune.jl` is the model for testing a post-processing stage (hand-build a net with `MM.add_vertex!`/`MM.add_element!`, transform, assert invariants). A new `test/test_smooth.jl` needs an `include` line added in `test/runtests.jl`.

## Decisions & conclusions
- Bounded change: one new exported function (e.g. `smooth`/`smooth_network` in a new `src/smooth.jl`), same-type-in/same-type-out so the Makie ext and any export path work unchanged. Subdivision must push one radius per new sub-element and one `Interior` tag + zero bcvalue per new vertex; endpoints keep their original tags/bcs.
- **Kyle's requirement (2026-08-25): first-derivative (C¹) continuity at shared vertices.** This rules out independent per-segment wiggling (kinks where segments meet). Decided mechanism: **per-vertex tangents + one cubic Hermite curve per element**, sampled at spacing `h` into sub-elements. Tangent assignment by vertex degree: degree 2 → one shared Catmull–Rom tangent `(B − A)/|B − A|` from the two neighbors (C¹ through the vertex by construction); degree 1 (inlets/outlets) → own chord direction; degree ≥ 3 (junctions) → full C¹ is geometrically impossible, so pair the two most collinear (or largest-radius) incident elements as a "through vessel" sharing a tangent (C¹ through the bifurcation), remaining branches attach at an angle (anatomically correct).
- **Where curvature comes from**: a Hermite cubic with chord-aligned end tangents is exactly straight, and the extruded capillary lattice has essentially no degree-2 interior nodes — every capillary node is a junction of vertical + horizontal elements (`src/extrude.jl:63-74`), so chains exist only in the venule trunk (`src/trees.jl:49-67`) and arteriole trees (`src/trees.jl:110-134`). Organic curvature therefore comes from **small random tilts of the end tangents** (seeded RNG) at junction/endpoint ends only; shared degree-2 tangents are computed, never randomized, preserving C¹ along vessels. Vertex positions never move → topology, junctions, and BC tags untouched. Knobs: max tilt angle, sample spacing `h`, `rng`.
- Junctions (degree ≥ 3) and tagged boundary vertices must stay fixed so connectivity and BCs are preserved. No adjacency/degree utility exists (Graphs.jl is NOT a dep — deps are DelaunayTriangulation, LinearAlgebra, Random, StaticArrays); build degree counts with a simple pass over `elements`. Union-find `_find!` in `src/prune.jl:4-10` is the only existing graph helper.
- Follow the element-type parametrization (commit `083a6a7`): `T<:Real` (not AbstractFloat, deliberately autodiff-capable), wrap literals as `T(x)`, and use the `_randn(rng, T)` idiom from `src/voronoi.jl:11` for any noise. Take an `rng::AbstractRNG` first argument like the other stages.

## What's left / next steps
1. Resolve with Kyle: (a) default values for the max tilt angle and sample spacing `h`; (b) whether `generate_network` calls smoothing automatically (config-gated) or it stays a separate user-invoked step; (c) through-vessel pairing rule at junctions — most-collinear vs largest-radius; (d) whether tilt should scale by vessel size (radius threshold is the only available signal). The mechanism itself (per-vertex tangents + Hermite, C¹ at degree-2 vertices, random end-tangent tilts) is settled — see Decisions.
2. Present the short design in chat (approach, files touched, tests) and get approval.
3. Implement `src/smooth.jl` (+ include/export in `src/MicrovascularMeshes.jl`), reusing `add_vertex!`/`add_element!`; keep it pure/out-of-place like prune.
4. Write `test/test_smooth.jl` modeled on `test/test_prune.jl` (invariant checks, endpoint/tag preservation, connectivity via the `ncomponents` helper, both TEST_TYPES) and add its include line; verbose test output per Kyle's rules; do NOT run tests unless told.
5. Optionally update an example script to show a smoothed render.

## Gotchas / constraints
- Radii are per **element**, tags/bcs per **vertex** — subdividing must maintain both invariants or plotting and any solver downstream silently breaks.
- The four generation stages are order-dependent; smoothing must run **after** `prune_disconnected` (it renumbers vertices).
- Plotting lives in the Makie extension only; smoothing stays in `src/`, Makie-free. The ext hardcodes `Vector{Float64}` for the color channel (`ext/…MakieExt.jl:22`) — harmless, but don't propagate that pattern.
- `voronoi_network` drops to Float64 internally for DelaunayTriangulation — irrelevant to smoothing, but don't copy that pattern; smoothing should stay generic in `T`.
- Kyle's testing rules: tests for every new function, verbose output, never run tests unless explicitly told.
- No secrets involved in this task.
