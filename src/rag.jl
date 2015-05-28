# Return a regional adjacency graph (of type Graph, as defined in 
# LightGraphs.jl). This is an undirected planar graph. The nodes correspond
# to segment labels. Edges join adjacent nodes only.
function adjacency_graph(labels::AbstractArray,
                         nlabels::Integer,
                         lab_means::Dict,
                         dep_means::Dict;
                         depth_weight::FloatingPoint=0.5)

    @assert 0.0 <= depth_weight <= 1.0

    g = Graph(nlabels)
    dep_dists = spzeros(nlabels, nlabels)
    lab_dists = spzeros(nlabels, nlabels)
    nr, nc = size(labels)
    borders = zeros(nr, nc)

    # Relative neighbor indices
    i4 = [-1,  0,  1,  0]
    j4 = [ 0, -1,  0,  1]

    for i = 1:nr
        for j = 1:nc
            for k = 1:4
                ni, nj = i+i4[k], j+j4[k]
                (1 <= ni <= nr) && (1 <= nj <= nc) || continue
                labels[i,j] != labels[ni, nj] || continue
                borders[i,j] = 1
                !has_edge(g, labels[i,j], labels[ni, nj]) || continue
                a, b = labels[i,j], labels[ni, nj]
                ca, cb = lab_means[a], lab_means[b]
                            
                rgb = convert(Color.RGB, ca)

                add_edge!(g, a, b)
                # Assign edge weights based on color and depth
                # Color distance
                cdiff = [ca.l - cb.l, ca.a - cb.a, ca.b - cb.b]
                lab_dists[a,b] = dot(cdiff, cdiff)
                # Depth distance
                dep_dists[a,b] = abs(dep_means[a] - dep_means[b])
            end
        end
    end

    lab_dists -= minimum(lab_dists)
    dep_dists -= minimum(dep_dists)
    lab_dists /= maximum(lab_dists)
    dep_dists /= maximum(dep_dists)

    g, ((1-depth_weight)*lab_dists + depth_weight*dep_dists), borders
    # g, 0.25*lab_dists + 0.75*dep_dists, borders
end

function cut_graph(g::Graph, edgewts::AbstractArray, thresh::FloatingPoint)
    g2 = Graph(nv(g))
    for e in edges(g)
        a, b = src(e), dst(e)
        wt = edgewts[a,b]
        if wt < thresh 
            add_edge!(g2, e)
        end
    end
    g2
end

function cut_graph!(g::Graph, edgewts::AbstractArray, thresh::FloatingPoint)
    for e in edges(g)
        a, b = src(e), dst(e)
        wt = edgewts[a,b]
        if wt > thresh 
            rem_edge!(g, e)
        end
    end
end

