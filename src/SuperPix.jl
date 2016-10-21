module SuperPix

using Images
using Colors
using LightGraphs

include("grids.jl")
include("means.jl")
include("prep.jl")
include("rag.jl")
include("recombine.jl")
include("slic.jl")
include("vis.jl")

export
    hexgrid,
    adjusted_grid,
    gradmag2,
    generalized_mean,
    color_means,
    color_moments,
    cluster_centroids,
    unpack_drgb_array,
    lab_array,
    preprocess_depth_array!,
    adjacency_graph,
    merged_superpixels,
    update_distances!,
    update_centroids!,
    reconnect,
    slic,
    graph_image,
    segment_borders
end
