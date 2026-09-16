module MicrovascularMeshesMakieExt

using MicrovascularMeshes
import Makie

function MicrovascularMeshes.plotnetwork(net::VascularNetwork; kwargs...)
    fig = Makie.Figure()
    ax = Makie.Axis3(fig[1, 1]; aspect = :data,
        title = "Idealized myocardial microvascular network",
        xlabel = "x", ylabel = "y", zlabel = "z")
    plotnetwork!(ax, net; kwargs...)
    return fig
end

function MicrovascularMeshes.plotnetwork!(ax, net::VascularNetwork;
        colormap = :viridis, linewidth = 1.5, show_nodes = true)
    verts = net.vertices
    nelem = length(net.elements)

    # Two endpoints per element; each segment is colored by its vessel radius.
    segments = Vector{Makie.Point3f}(undef, 2 * nelem)
    segment_radius = Vector{Float64}(undef, 2 * nelem)
    for (k, element) in enumerate(net.elements)
        segments[2k-1] = Makie.Point3f(verts[element[1]]...)
        segments[2k] = Makie.Point3f(verts[element[2]]...)
        segment_radius[2k-1] = net.radii[k]
        segment_radius[2k] = net.radii[k]
    end
    Makie.linesegments!(ax, segments;
        color = segment_radius, colormap = colormap, linewidth = linewidth)

    if show_nodes
        for (tag, color, markersize) in (
            (NoFlow, :gray, 4),
            (Outlet, :dodgerblue, 14),
            (Inlet, :crimson, 14),
            (InletVenule, :orange, 14),
        )
            idx = findall(==(tag), net.node_tags)
            isempty(idx) && continue
            points = [Makie.Point3f(verts[i]...) for i in idx]
            Makie.scatter!(ax, points; color = color, markersize = markersize)
        end
    end
    return ax
end

end # module MicrovascularMeshesMakieExt
