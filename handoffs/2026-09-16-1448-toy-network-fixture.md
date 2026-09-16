---
slug: toy-network-fixture
created: 2026-09-16-1448
status: done  # implemented as toy_network in commit b776ad0 (radial fans spaced π/copies, not 2π/copies); :stacked later dropped and :volumetric made the default
---

# Handoff: `toy_network` — hand-built arteriole → capillary bed → venule fixture

## Goal / why this matters
Add one exported function, `toy_network`, that hand-builds a small deterministic `VascularNetwork` shaped like the textbook diagram Kyle supplied: a bifurcating arterial tree fanning out from a single inlet on the left, an anastomotic capillary bed in the middle, and a mirror-image venous tree converging to a single outlet on the right. `generate_network` is stochastic and produces thousands of elements, so it is a poor fixture for testing solvers, exporters, or `smooth_network`; this gives tests and demos a tiny network whose counts, radii, and symmetries are known in closed form.

## Background & current state
Brainstorming and plan mode completed on 2026-09-16; **design is settled and Kyle answered the three requirement questions** (below), but he parked the work before approving the implementation plan. No code was written. Before implementing, present the short design (this doc's "Decisions" section) in chat and get an explicit yes — the brainstorming approval gate was never passed.

Kyle's decisions:
1. **All three 3D readings are wanted, plus flat** — he is not sure which 3D shape he wants, so build `:flat`, `:volumetric`, `:stacked`, and `:radial` and let him look at a 4-panel figure. Expect him to delete variants afterwards; keep them cheap to remove.
2. **Lives in `src/`, exported** (like `Graphs.smallgraph`), so tests, `examples/`, docs, and downstream solver packages can all reach it.
3. **Capillary bed is an anastomotic ladder** (arterial leaf k → venous leaf k, plus arterial leaf k → venous leaf k+1), producing loops; not one-to-one pairing, no keyword to switch it off.

Open naming choice: `toy_network` (proposed) vs `example_network`. Ask when presenting the design.

## Key files / locations
- `src/types.jl:92-102` — `VascularNetwork{T<:Real}`: `vertices::Vector{SVector{3,T}}`, `elements::Vector{SVector{2,Int}}`, `radii::Vector{T}` (**per element**), `node_tags::Vector{NodeTag}` and `bcvalues::Vector{T}` (**per vertex**). No inner constructor; invariants are documented only. `T<:Real` deliberately (autodiff), wrap literals as `T(x)`.
- `src/types.jl:15-21` — `NodeTag` enum: `Interior`, `NoFlow`, `Inlet`, `Outlet`, `InletVenule`. There is no arterial/venous/capillary label; radius is the only vessel-class signal.
- `src/types.jl:116` `add_vertex!` (tags `Interior`, bc 0), `src/types.jl:129` `add_element!` (pushes element + radius), `src/trees.jl:4` `set_bc!`. Build with these; none are exported (tests use `const IM = IdealizedMicrovasculature`).
- BC conventions in `src/trees.jl:59-60,147`: `Inlet → 1e-3` (velocity), `Outlet → 2000` (pressure). Package trunks are a `start → finish` stub so tagged vertices have degree 1.
- `src/IdealizedMicrovasculature.jl:8-21` — export lines and `include` order (domains, types, voronoi, extrude, trees, prune, smooth, generate). Add `include("toy.jl")` after `generate.jl` and `toy_network` to the `generate_network, smooth_network` export line.
- `src/smooth.jl` — `smooth_network(rng, net; max_tilt, spacing)`; the fixture's curvature comes from composing with this, not from the builder.
- `ext/IdealizedMicrovasculatureMakieExt.jl:15-46` — `plotnetwork!` reads `vertices`, `elements`, `radii` (color channel), `node_tags` (scatter for `Inlet`/`Outlet`/`InletVenule`/`NoFlow`). Positive radii + tags are all a plottable fixture needs. Uses `Axis3(aspect = :data)`, which degenerates on a zero z-extent (the flat layout) — give the flat panel an explicit `aspect` and a top-down `elevation` in the demo.
- `test/runtests.jl` — one outer `@testset`, one `include` per module, `TEST_TYPES = (Float64, Float32)`, and the `ncomponents(net)` union-find helper. `test/test_smooth.jl` is the model: tiny hand-built nets, exact counts, invariant lengths, tag preservation, `@test_throws ArgumentError` validation, both element types.
- `examples/demo_smoothing.jl` — self-bootstrapping demo script pattern (temp env, `Pkg.develop` the repo, `CairoMakie`, save PNG next to the script). Copy its shape for `examples/demo_toy.jl`.
- `docs/pages/api.qmd` iterates `names(IdealizedMicrovasculature)` and splices docstrings, so the new export needs a house-style docstring (`### Keyword Arguments`, `### Examples`) and nothing else for docs.
- No Graphs.jl dependency; do not add one. Deps are DelaunayTriangulation, LinearAlgebra, Random, StaticArrays.

## Decisions & conclusions

**API — one function, one layout keyword.** All four variants share one tree builder and differ only in leaf placement and how copies join, so a `layout::Symbol` is enough; four layout structs would be more surface than warranted for variants Kyle may delete.

```julia
toy_network(::Type{T}=Float64;
    layout=:flat,          # :flat | :volumetric | :stacked | :radial
    generations=3,         # bifurcations per tree; 2^generations leaves
    extent=200,            # inlet-to-outlet distance along x
    width=100,             # transverse spread of the leaves
    capillary_radius=2.6,  # matches NetworkConfig.radius_default_small_vessels
    layers=3,              # :stacked only
    copies=4,              # :radial only
) -> VascularNetwork{T}
```
Type parameter follows the `zeros(Float32, n)` idiom. `extent` rather than `length` so `Base.length` is not shadowed inside the body. Pure and deterministic — no rng.

**Geometry along x** (fixed fractions of `extent`, not exposed): trunk 0.15, arterial fan 0.25, capillary bed 0.2, venous fan 0.25, trunk 0.15. Inlet vertex at x = 0, arterial hub at 0.15·extent, generation g node at hub + 0.25·extent·g/G, bed occupies the middle 0.2·extent, venous side is the exact mirror x ↦ extent − x.

**Tree as a binary path — the mechanism that unifies the layouts.** Leaf k (k−1 written in G bits, most-significant bit = generation 1) sits at transverse offset Σ_g s_g·δ_g·â(g), where s_g = ±1 is the bit for generation g and â(g) is that generation's split axis. A node at generation g is the partial sum of its path, so every internal node is the centroid of its descendant leaves and the tree fans out as a lens. δ_g halves each time an axis is reused: δ_g = (width/2) / 2^(splits so far on that axis, counting this one). Under this one rule:
- `:flat` — â(g) = ŷ always. Leaves evenly spaced in y at pitch width/2^G, all z = 0. This is the picture.
- `:volumetric` — â(g) alternates ŷ, ẑ, ŷ, … Leaves form a 2^⌈G/2⌉ × 2^⌊G/2⌋ grid in y–z. Consecutive leaf indices follow a Z-order curve, so a few ladder capillaries are short diagonals; document it, don't fight it.
- `:stacked` — the flat network copied `layers` times at z = (ℓ−1)·(width/2^G) (the leaf pitch, so vertical links are capillary-scale). Corresponding arterial and venous **leaves** in adjacent layers are joined by capillary-radius vertical elements, mirroring `extrude_capillary!`. Each layer keeps its own `Inlet`/`Outlet` (the generated network is multi-inlet too).
- `:radial` — one inlet and one outlet on the x axis; the flat pattern is rotated about x into `copies` azimuthal copies (ŷ ↦ ŷ cos φ + ẑ sin φ, φ = 2π(c−1)/copies) sharing trunk and hub. Hub degree is 2·copies + 1.

**Capillary bed (anastomotic ladder).** n = 2^G leaves per tree: arterial leaf k → venous leaf k for all k, plus arterial leaf k → venous leaf k+1 for k < n; 2n − 1 capillaries at `capillary_radius`. Loops close through the trees, as in the picture. Per layer in `:stacked`, per copy in `:radial`.

**Radii by Murray's law** (r_parent³ = Σ r_child³), driven upward from `capillary_radius`: the vessel entering a generation-g node carries r_g = capillary_radius·2^((G+1−g)/3); the trunk carries r_0 (× copies^(1/3) for `:radial`). One knob, physiological, and it makes the big-vessel/capillary distinction visible in the radius colormap. Murray holds exactly at every tree junction; the two ladder ends (arterial leaf n, venous leaf 1) have one capillary instead of two, which is acceptable.

**Tags/BCs** follow `src/trees.jl`: arterial inlet `Inlet`/1e-3, venous outlet `Outlet`/2000, everything else `Interior`/0. Inlet and outlet vertices are degree 1 (that is what the trunk stub is for). No `NoFlow` tags.

**Closed-form counts** for tests (n = 2^G): flat vertices = 2^(G+2), flat elements = 3·2^(G+1) − 3 (check G=1: 8 vertices, 9 elements). Stacked = layers × flat + (layers−1)·2n vertical links. Radial vertices per side = 2 + copies·(2^(G+1) − 2), elements per side = 1 + copies·(2^(G+1) − 2), capillaries = copies·(2n − 1).

## What's left / next steps
1. Present the design above in chat (approach, files touched, tests); confirm the name (`toy_network` vs `example_network`); get an explicit yes.
2. Write `src/toy.jl`: `toy_network` + private helpers (`_tree_positions`, `_add_tree!`, `_add_ladder!`, per-layout assembly). Validate `generations ≥ 1`, `layers ≥ 1`, `copies ≥ 1`, positive sizes, known `layout`; `ArgumentError` otherwise. House-style docstring showing composition with `smooth_network`.
3. Wire `include`/`export` in `src/IdealizedMicrovasculature.jl`; add a Key-files line for `src/toy.jl` to `CLAUDE.md`.
4. Write `test/test_toy.jl` (+ `include("test_toy.jl")` after `test_smooth.jl` in `test/runtests.jl`), every layout × `TEST_TYPES`: length invariants; exact counts from the closed forms for G = 1, 2, 3; one `Inlet`/`Outlet` per layer (one total for `:radial`) with BC values 1e-3/2000 and degree 1; all other tags `Interior`/0; `ncomponents == 1`, indices in range, no self-loops or duplicates; Murray's law at every tree junction; exactly the expected number of `capillary_radius` elements; `:flat` all z == 0 with evenly spaced leaves; `:volumetric` spread in both y and z; `:stacked` exactly `layers` distinct z; `:radial` invariant under rotation by 2π/copies about x; every layout invariant under x ↦ extent − x; `Float32` preserved; `smooth_network(MersenneTwister(0), toy_network())` keeps one component and the tagged vertices; `@test_throws ArgumentError` for bad inputs. Verbose testset names carrying `layout` and `T`. **Write tests, do not run them unless Kyle says so.**
5. Write `examples/demo_toy.jl` (2×2 panels of the four layouts via `plotnetwork!`, flat panel top-down with explicit `aspect`) and render `examples/demo_toy.png`; this figure is how Kyle will pick which 3D variant(s) to keep.
6. Verify in the Julia daemon: build each layout for G = 1..3 and compare counts to the closed forms; smooth the flat network and confirm junctions/tags unchanged.

## Gotchas / constraints
- Radii are per element, tags/bcs per vertex — every `add_element!` needs a radius, every `add_vertex!` gets `Interior`/0 by default; only the two roots per tree are retagged with `set_bc!`.
- Keep `src/` Makie-free; plotting is the extension's job. Do not hardcode `Float64` anywhere (the ext does for its color buffer; don't copy that).
- Follow the loaded style rules: `T(x)` literals, `map`/`foreach` over comprehensions where natural, no spaces around `=` in kwargs, unicode `≤`/`≥`/`∈`, indented one-per-line argument blocks when a signature breaks.
- Kyle's rules: tests for every new function, verbose test output, never run tests unless told; commit messages `<type>(<scope>): <subject>` with the `agent-footer.sh --trailers` trailers; open any markdown deliverable in Peekaboo afterwards (`peekaboo "<path>"`).
- The plan-mode file at `~/.claude/plans/lets-add-a-helper-vivid-rabbit.md` holds the same design; this handoff supersedes it.
- No secrets involved.
