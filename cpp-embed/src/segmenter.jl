include("/Users/adare/repos/superpix/slic.jl")

function unpack_drgb_array(pbuff::Array{UInt32, 1}, nrows::Int64, ncols::Int64)
    
    pbuff_2d = reshape(pbuff, ncols, nrows)

    # Unpack UInt32s into LAB color and depth arrays.
    dep  = Array(Float32, ncols, nrows)
    alab = Array(Color.Lab{Float32}, ncols, nrows)

    for i = 1:ncols
        for j = 1:nrows
            p = UInt32(pbuff_2d[i,j])
            b = ((p >>  0) & 0xFF)/255
            g = ((p >>  8) & 0xFF)/255
            r = ((p >> 16) & 0xFF)/255
            d = ((p >> 24) & 0xFF)/255
            dep[i,j] = d
            alab[i,j] = convert(Color.Lab{Float32}, Color.RGB(r, g, b))
        end
    end
    dep, alab
end

function preprocess_depth_array!(dep::AbstractArray)
    mindep = minimum(dep[dep .> 0])
    dep = clamp(dep, mindep, maximum(dep)) - mindep
    dep /= maximum(dep)
    # println("min, max = $(minimum(dep)), $(maximum(dep))")

    # This works:
    # depth_img = grayim(dep)
    # imwrite(depth_img, "depth.jpg")
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

function segment_drgb(pbuff::Array{UInt32, 1}, nrows::Integer, ncols::Integer)

    dep, alab = unpack_drgb_array(pbuff, nrows, ncols)
    preprocess_depth_array!(dep)

    imlab = Image(alab, IMcs="sRGB", spatialorder=["x","y"], pixelspacing=[1,1])

    k, m = 1000, 10
    @time labels, nlabels = slic(imlab, k, m)

    lab_means, color_superpix = color_means(imlab, labels, nlabels)
    dep_means, depth_superpix = color_means(dep, labels, nlabels)

    nr, nc = size(labels)
    graph, edgewts, borders = adjacency_graph(labels, nlabels, lab_means, dep_means)
    centroids = cluster_centroids(labels, nlabels)

    println("Adjacency graph has $(nv(graph)) vertices and $(ne(graph)) edges.")
    cut_graph!(graph, edgewts, 0.02)

    seg_labels, seg_borders = merged_superpixels(labels, nlabels, graph)

    centroid_img = zeros(labels)
    for c in 1:nlabels
        row, col = centroids[c, 1], centroids[c, 2]
        if (1 <= row <= nr) && (1 <= col <= nc)
            centroid_img[row, col] = 1
        end
    end

    graph_edges = graph_image(graph, centroids, nr, nc)

    graphcuts = Overlay((borders', graph_edges', centroid_img'),
                        (Color.RGB(0.2,0.2,0.2), Color.RGB(0,0.5,1), Color.RGB(1,1,0)),
                        ((0,1), (0,1), (0,1))
                        )

    imwrite(grayim(seg_borders), "seg_borders.jpg")
    imwrite(graphcuts, "graphcuts.jpg")
    # imwrite(segs, "segs.jpg")
end

function test()
    pbuff = readdlm("../pbuff.txt", UInt32)
    println(size(pbuff))
    nrows, ncols = 290, 700
    # nrows, ncols = 171, 361

    segment_drgb(pbuff[:], nrows, ncols)
end
