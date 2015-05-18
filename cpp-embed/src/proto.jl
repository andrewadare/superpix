using Images

include("/Users/adare/repos/superpix/slic.jl")

function main()
    pbuff = readdlm("../pbuff.txt", UInt32)
    nrows, ncols = 171, 361
    pbuff = reshape(pbuff, ncols, nrows)

    # Unpack UInt32s into LAB color and depth arrays.
    dep  = Array(Float32, ncols, nrows)
    alab = Array(Color.Lab{Float32}, ncols, nrows)

    for i = 1:ncols
        for j = 1:nrows
            p = UInt32(pbuff[i,j])
            b = ((p >>  0) & 0xFF)/255
            g = ((p >>  8) & 0xFF)/255
            r = ((p >> 16) & 0xFF)/255
            d = ((p >> 24) & 0xFF)/255

            dep[i,j] = d
            alab[i,j] = convert(Color.Lab{Float32}, Color.RGB(r, g, b))
        end
    end

    # This works:
    # depth_img = grayim(dep)
    # imwrite(depth_img, "depth.jpg")

    imlab = Image(alab)
    imlab["IMcs"] = "sRGB"
    imlab["spatialorder"] = ["x","y"]
    imlab["pixelspacing"] = [1,1]
    show(imlab)

    k, m = 1000, 10
    @time labels, nlabels = slic(imlab, k, m)

    lab_means, color_superpix = color_means(imlab, labels, nlabels)
    dep_means, depth_superpix = color_means(dep, labels, nlabels)

    nr, nc = size(labels)
    graph, edgewts, borders = adjacency_graph(labels, nlabels, lab_means, dep_means)
    centroids = cluster_centroids(labels, nlabels)

    println("Adjacency graph has $(nv(graph)) vertices and $(ne(graph)) edges.")
    for e in edges(graph)
        a, b = src(e), dst(e)
        wt = edgewts[a,b]
        thresh = 0.18
        flag = ""
        if wt > thresh 
            flag = "<-- remove"
            rem_edge!(graph, e)
        end
    end

    centroid_img = zeros(labels)
    for c in 1:nlabels
        row, col = centroids[c, 1], centroids[c, 2]
        if (1 <= row <= nr) && (1 <= col <= nc)
            centroid_img[row, col] = 1
        end
    end

    graph_edges = graph_image(graph, centroids, nr, nc)

    segs = Overlay((borders', graph_edges', centroid_img'), 
                   (Color.RGB(0.2,0.2,0.2), Color.RGB(0,0.5,1), Color.RGB(1,1,0)),
                   ((0,1), (0,1), (0,1))
                   )


    imwrite(segs, "segs.jpg")
    # imwrite(imlab, "imlab.jpg")
end

main()
