using Images
using Color
using LightGraphs
using Reactive, Interact

include("grids.jl")
include("means.jl")
include("io.jl")
include("rag.jl")
include("recombine.jl")
include("slic.jl")
include("vis.jl")
include("gradients.jl")


# Run over a pair of images named e.g. "../input/color_1.jpg" and 
# "../input/depth_1.png". The results will be written to ../output.
#
# Example usage in julia REPL:
# julia> include("color_depth_seg.jl"); @time for i = 100:101 @time seg(i) end

function seg(idx::Integer)
    println("Segmenting image $idx")

    row_range = 91:380      # Crop input images to this region
    col_range = 151:750
    k = 1000                # Number of requested superpixels
    m = 10                  # Cluster compactness parameter
    w = 0.5                 # Depth vs. color weighting for RAG edges
    t = 0.07                # Graph cut threshold

    imlab = lab_image("../input/color_$idx.jpg", row_range, col_range)
    dep = depth_image("../input/depth_$idx.png", row_range, col_range)
 
    depgrad = depth_image_grad(dep)
    labgrad = color_image_grad(imlab)

    nr, nc = size(imlab)
    seeds = adjusted_grid(hexgrid(nr, nc, k), depgrad+labgrad)

    # Do superpixel segmentation
    sp_labels, nsp = slic(imlab, seeds, k, m, niter=5)
    sp_centroids = cluster_centroids(sp_labels, nsp)

    # Mean color values over superpixels
    lab_means, color_superpix = color_means(imlab, sp_labels, nsp)
    dep_means, depth_superpix = color_means(dep,   sp_labels, nsp)

    # Compute a regional adjacency graph from superpixels, then cut it
    graph, edgewts, borders = 
    adjacency_graph(sp_labels, nsp, lab_means, dep_means, depth_weight=w)
    graph = cut_graph(graph, edgewts, t)

    # Make an image of the graph including superpixel centroids
    graph_edges = graph_image(graph, sp_centroids, nr, nc)

    # Combine superpixels according to the cut adjacency graph
    seg_labels, nsegments = merged_superpixels(sp_labels, graph)

    # Make an image of colored segment patches
    seg_lab_means, color_segments = color_means(imlab, seg_labels, nsegments)
    
    # Overlay segment boundaries on (a copy of) the input color image
    seg_img = segment_overlay(seg_labels, imlab, color_segments,
                              cluster_centroids(seg_labels, nsegments))

    imwrite(convert(Image{Color.RGB}, imlab),          "../output/color_in_$idx.jpg")
    imwrite(grayim(dep),                               "../output/depth_in_$idx.jpg")
    imwrite(convert(Image{Color.RGB}, seg_img),        "../output/overlay_$idx.jpg")
    imwrite(convert(Image{Color.RGB}, color_segments), "../output/color_segments_$idx.jpg")
    imwrite(convert(Image{Color.RGB}, color_superpix), "../output/color_superpix_$idx.jpg")
    imwrite(grayim(depth_superpix),                    "../output/depth_superpix_$idx.jpg")
    imwrite(grayim(graph_edges),                       "../output/graph_edges_$idx.jpg")

    # grayim(depgrad)
    # grayim(gradmean)
end
