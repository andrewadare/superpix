# Return a regional adjacency graph (of type Graph, as defined in
# LightGraphs.jl). This is an undirected planar graph. The nodes correspond
# to segment labels. Edges join adjacent nodes only.
function adjacency_graph(labels::AbstractArray,
                         nlabels::Integer,
                         lab_means::Dict,
                         dep_means::Dict;
                         depth_weight::AbstractFloat=0.5)

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

                rgb = convert(Colors.RGB, ca)

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

# Compute the integrated intensity of image gradient magnitude between nodes
# and return in a sparse matrix.
function grad_weights(g::Graph, centroids::AbstractArray, grad_img::AbstractArray)
    nr, nc = size(grad_img)
    wts = spzeros(Float64, nv(g), nv(g))
    for e in edges(g)
        a, b = src(e), dst(e)
        r1, c1 = centroids[a, 1], centroids[a, 2]
        r2, c2 = centroids[b, 1], centroids[b, 2]
        if (1 <= r1 <= nr) && (1 <= c1 <= nc) && (1 <= r2 <= nr) && (1 <= c2 <= nc)
            dist = hypot(r2-r1, c2-c1)
            theta = atan2(r2-r1, c2-c1)
            i, d = 0, 0
            while d < dist
                dr = i*sin(theta)
                dc = i*cos(theta)
                wts[a,b] += grad_img[round(Int, r1+dr), round(Int, c1+dc)]
                d = hypot(dr, dc)
                i += 1
            end
        end
    end
    wts
end

function cut_graph(g::Graph, edgewts::AbstractArray, thresh::AbstractFloat)
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

function cut_graph!(g::Graph, edgewts::AbstractArray, thresh::AbstractFloat)
    for e in edges(g)
        a, b = src(e), dst(e)
        wt = edgewts[a,b]
        if wt > thresh
            rem_edge!(g, e)
        end
    end
end

