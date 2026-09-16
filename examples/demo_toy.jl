#!/usr/bin/env julia
#
# Standalone demo of toy_network: the three layouts of the hand-built
# arteriole → capillary ladder → venule fixture, side by side.
#
# Run it directly from a fresh checkout — no environment setup required:
#
#     julia examples/demo_toy.jl
#
# The script bootstraps its own throwaway environment, develops the
# MicrovascularMeshes package from this repository, and adds CairoMakie
# (a headless, no-GPU backend that writes a static PNG). Swap CairoMakie for
# GLMakie below for an interactive, rotatable window. Packages are cached
# globally, so only the first run is slow.

using Pkg

repo_root = dirname(@__DIR__)
Pkg.activate(; temp=true)
Pkg.develop(; path=repo_root)
Pkg.add("CairoMakie")

using CairoMakie
using MicrovascularMeshes

#--------------------------------------------------------------------------------# Build

# layout, panel title, and grid position: the wide flat picture spans the top row
layouts = (
    (:flat, "flat — the textbook picture (top-down)", (1, 1:2)),
    (:volumetric, "volumetric (default) — splits alternate y and z", (2, 1)),
    (:radial, "radial — 4 flat fans rotated about x", (2, 2)),
)
networks = map(((layout, _, _),) -> toy_network(; layout), layouts)

foreach(layouts, networks) do (layout, _, _), net
    println(rpad(":$layout", 12), "$(nvertices(net)) vertices, $(nelements(net)) elements")
end

#--------------------------------------------------------------------------------# Visualize

# Axis3(aspect=:data) collapses on a zero-extent axis (z ≡ 0 for the flat layout), so the
# flat panel is drawn top-down in a 2D axis and the others get a box sized from the data.
function box_aspect(net)
    lo = reduce((a, b) -> min.(a, b), net.vertices)
    hi = reduce((a, b) -> max.(a, b), net.vertices)
    spans = max.(hi - lo, 0.05 * maximum(hi - lo))
    return Tuple(spans ./ maximum(spans))
end

function panel!(fig, pos, net, title; top_down=false)
    ax = if top_down
        Axis(fig[pos...]; aspect=DataAspect(), title, xlabel="x", ylabel="y")
    else
        Axis3(fig[pos...]; aspect=box_aspect(net), title, azimuth=1.3π, elevation=π / 7,
            xlabel="x", ylabel="y", zlabel="z")
    end
    plotnetwork!(ax, net; linewidth=3)
    return ax
end

fig = Figure(size=(1500, 1300))
Label(fig[0, 1:2], "toy_network layouts — color is vessel radius (Murray's law), " *
    "red = inlet, blue = outlet"; fontsize=22, font=:bold)

foreach(layouts, networks) do (layout, title, pos), net
    panel!(fig, pos, net, title; top_down=layout ≡ :flat)
end

snapshot = joinpath(@__DIR__, "demo_toy.png")
save(snapshot, fig; px_per_unit=2)
println("\nSaved visualization to $snapshot")
