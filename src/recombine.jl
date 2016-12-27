function merged_superpixels(labels::AbstractArray, graph::Graph)
    nr, nc = size(labels)

    @assert nr > 0 && nc > 0 "nr, nc = $nr, $nc"

    seglabels = zeros(Int, nr, nc)

    g = Graph(nv(graph))
    for v in vertices(graph)
        if degree(g, v) > 0
            add_vertex!(g, v)
        end
    end
    for e in edges(graph)
        if has_vertex(g, src(e)) && has_vertex(g, dst(e))
            add_edge!(g,e)
        end
    end

    label = 0
    for s in vertices(g)
        vv = visited_vertices(g, DepthFirst(), s)

        if length(vv) > 2
            label += 1
            for v in vv
                indices = find(labels .== v)

                for i in indices
                    if seglabels[i] == 0
                        seglabels[i] = label
                    end
                end
            end
        end

        for v in vv
            println("***** ", v)
            v_neighbors = []
            try
                v_neighbors = neighbors(g, v)
                for neighbor in v_neighbors
                    # println("g.fadjlist[v], v, neighbor = $(g.fadjlist[v]), $v, $neighbor")
                    # if g.fadjlist
                    rem_edge!(g, v, neighbor)
                end
            catch
                println("g.fadjlist[v], v, v_neighbors = $(g.fadjlist[v]), $v, $v_neighbors")
            end
        end
    end

    seglabels, label
end
