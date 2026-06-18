#--------------------------------------------------------------------------------# Pruning

# Union-find root of `x` with path halving.
function _find!(parent::Vector{Int}, x::Int)
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    return x
end

"""
    prune_disconnected(net) -> VascularNetwork

Drop capillaries that are not part of the main network: return a copy of `net` keeping only
its largest connected component, with vertices, elements, radii, and boundary tags/values
re-indexed consistently. Vessels reachable from the inlet/outlet trees live in that
component, so the discarded pieces are exactly the unperfusable, free-floating fragments.

A fully connected network is returned unchanged (the identical object).
"""
function prune_disconnected(net::VascularNetwork{T}) where {T}
    n = nvertices(net)
    n == 0 && return net

    parent = collect(1:n)
    for e in net.elements
        a, b = _find!(parent, e[1]), _find!(parent, e[2])
        a != b && (parent[a] = b)
    end

    sizes = zeros(Int, n)
    for v in 1:n
        sizes[_find!(parent, v)] += 1
    end
    keep_root = argmax(sizes)

    keep = [_find!(parent, v) == keep_root for v in 1:n]
    all(keep) && return net

    newindex = zeros(Int, n)
    count = 0
    for v in 1:n
        keep[v] && (newindex[v] = (count += 1))
    end

    elements = SVector{2,Int}[]
    radii = T[]
    for (k, e) in enumerate(net.elements)
        keep[e[1]] || continue            # both endpoints share a component
        push!(elements, SVector{2,Int}(newindex[e[1]], newindex[e[2]]))
        push!(radii, net.radii[k])
    end
    return VascularNetwork{T}(net.vertices[keep], elements, radii,
                              net.node_tags[keep], net.bcvalues[keep])
end
