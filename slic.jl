using Images
# using Color

# Compute a hexagonal grid over a 2D region of size h x w.
# An n x 2 array of row,col positions is returned, where n <= k is the number
# of grid points.
# For a nice reference, see http://www.redblobgames.com/grids/hexagons/
function hexgrid(h::Integer, w::Integer, k::Integer)

    step = sqrt(h*w/k)
    xoff, yoff = round(Int, step/2), round(Int, step/2)
    grid = zeros(Int, k, 2)
    r = Int(0)
    n = Int(0)

    for i = 0:h
        row = round(Int, i*step + yoff)
        row > h && break
        for j = 0:w
            col = round(Int, j*step + xoff<<(r&1))
            col > w && break
            n += 1
            grid[n,:] = [row,col]
        end
        r += 1
    end

    grid[1:n,:]
end

# Compute squared magnitude of image gradient
function gradmag2(img::AbstractArray, method::String)
    gx, gy = imgradients(img, method)
    g = magnitude(gx, gy)

    sz = size_spatial(g)
    g2 = Array(Float64, sz)

    for i = 1:sz[1]
        for j = 1:sz[2]
            g2[i,j] = dot(g[:,i,j], g[:,i,j])
        end
    end

    # Return an array with the data of g2 and the properties of img.
    shareproperties(img, g2)
end

# Return grid points nudged to the minimum-gradient position in a 3x3 window
function adjusted_grid(grid::AbstractArray, grad_mag::AbstractArray)

    nr, nc = size(grad_mag)
    newgrid = copy(grid)

    for n = 1:size(grid,1)

        # Make a 3x3 window around the grid point i,j.
        i,j = grid[n,:]
        rows = max(1,i-1) : min(nr,i+1)
        cols = max(1,j-1) : min(nc,j+1)
        window = sub(grad_mag, rows, cols)

        # Find pixel position in window with the smallest gradient.    
        wi, wj = ind2sub(size(window), indmin(window))

        newgrid[n,:] = [rows[1] - 1 + wi, cols[1] - 1 + wj]
    end

    newgrid
end

function find_connected_labels(labels::AbstractArray, 
                               linear_label_index::Int,
                               segment::AbstractArray)
    npix = 1
    nr, nc = size(labels)
    i,j = ind2sub(size(labels), linear_label_index)

    # Relative neighbor indices
    i4 = [-1,  0,  1,  0]
    j4 = [ 0, -1,  0,  1]

    # Add this pixel to the list of connected pixels.
    segment[linear_label_index] = 1

    for k = 1:4
        ni, nj = i+i4[k], j+j4[k]
        if (1 <= ni <= nr) && (1 <= nj <= nc) && labels[ni,nj] == labels[i,j]
            neighbor_idx = sub2ind(size(labels), ni, nj)
            if segment[neighbor_idx] == 0
                npix += find_connected_labels(labels, neighbor_idx, segment)
            end
        end
    end
    npix
end

function reconnect(labels::AbstractArray, min_cluster_size::Integer)
    newlabels = zeros(labels)
    nr, nc = size(labels)

    # Relative neighbor indices
    i4 = [-1,  0,  1,  0]
    j4 = [ 0, -1,  0,  1]

    iseg = zeros(Int, nr*nc)
    jseg = zeros(Int, nr*nc)

    # # Linear indices of connected labels in one segment
    # segment = zeros(Int, nr*nc)
    # segsize = 0

    label = 1
    adjlabel = 1

    for i = 1:nr
        for j = 1:nc

            # Only enter the loop body if we are starting a new segment.
            newlabels[i,j] == 0 || continue
            newlabels[i,j] = label

            # Find an adjacent newlabel for later merging of small segments.
            for k = 1:4
                ni, nj = i+i4[k], j+j4[k]
                if (1 <= ni <= nr) && (1 <= nj <= nc)
                    if newlabels[ni, nj] > 0
                        adjlabel = newlabels[ni, nj]
                    end
                end
            end

            # segsize = 0
            # idx = newlabels[sub2ind(size(newlabels), i, j)]
            # println("idx: ", idx)
            # segsize = find_connected_labels(labels, idx, segment)
            # println("segsize: ", segsize)
            # for l=1:segsize
            #     println(segment[l])
            # #     newlabels[segment[l]] = label
            # end

            iseg[1], jseg[1] = i, j
            c, count = 0, 1
            while c < count
                c += 1
                for k = 1:4
                    ni, nj = iseg[c]+i4[k], jseg[c]+j4[k]
                    if (1 <= ni <= nr && 1 <= nj <= nc) && labels[i,j] == labels[ni,nj] && newlabels[ni,nj] == 0
                        count += 1
                        iseg[count], jseg[count] = ni, nj
                        newlabels[ni, nj] = label
                    end
                    # println("$ni, $nj: label $label newlabel $(newlabels[ni, nj])")
                end
            end

            # # Reassign small clusters to adjlabel & decrement label count.
            # if count < min_cluster_size
            #     for c = 1:count
            #         # println("iseg[$c],jseg[$c] = $(iseg[c]), $(jseg[c])")
            #         newlabels[iseg[c],jseg[c]] = adjlabel
            #     end
            #     label -= 1
            # end

            label += 1
        end
    end

    newlabels, label
end

# TODO(OPT) split image and use @parallel
function slic(img::AbstractArray, k::Integer, m::Integer)

    tic()
    nr, nc = size(img)

    # Seed the clusters
    grid = hexgrid(nr, nc, k)

    # Compute (squared) image gradient magnitude
    grad2 = gradmag2(img, "prewitt")

    # Shift grid points slightly away from edges
    centers = adjusted_grid(grid, grad2)
    nclusters = size(centers, 1)

    # Cluster labels for each pixel
    labels = zeros(Int, nr, nc)
    
    # Half-size of clustering window
    s = round(Int, sqrt(nr*nc/k)) - 2
    area = s*s

    max_npix = (2s+2)^2
    cluster_rows = Array(Int, max_npix, nclusters)
    cluster_cols = Array(Int, max_npix, nclusters)
    npix = zeros(Int, nclusters)
    m2 = m*m

    # Squared pixel-to-centroid Euclidean distance in LAB color space
    # color_d2 = zeros(Float64, nr, nc)


    # maxlab = 10*10*ones(Float64, nclusters) # Adaptive m-like parameter

    print("Before loop...")
    toc()

    for iter = 1:1

        # Distances from pixels to cluster centers, initialized to ∞
        d = convert(Float64, Inf)*ones(Float64, nr, nc)

        ########### TODO: put this in a update_distances!() function ###########
        tic()
        for n = 1:nclusters

            # Centroid of cluster n
            i,j = centers[n,:]

            # Clustering window (2s x 2s) centered at point i,j for seed n
            rows = max(1, i-s) : min(nr, i+s)
            cols = max(1, j-s) : min(nc, j+s)

            # Loop through clustering window
            for r in rows
                for c in cols
                    
                    # Squared pixel-to-centroid Euclidean distance
                    sdiff = [r,c] - [i,j]
                    pixel_d2 = dot(sdiff, sdiff)

                    # Squared color distance in LAB color space
                    # pc and cc are the pixel and centroid ColorValues
                    pc, cc = img[r,c], img[i,j]
                    cdiff = [pc.l - cc.l, pc.a - cc.a, pc.b - cc.b]
                    color_d2 = dot(cdiff, cdiff)
                    
                    dist = pixel_d2/area + color_d2/m2

                    # Assign new distances and labels
                    if dist < d[r,c]
                        d[r,c] = dist
                        labels[r,c] = n

                        # Store pixel position
                        npix[n] += 1
                        cluster_rows[npix[n], n] = r
                        cluster_cols[npix[n], n] = c
                    end

                end
            end
        end
        print("Dist calcs $iter...")
        toc()
        ########### TODO: put this in a update_distances!() function ###########

        ########### TODO: put this in a update_centroids!() function ###########
        tic()
        for n = 1:nclusters
            rows = cluster_rows[1:npix[n], n]
            cols = cluster_cols[1:npix[n], n]
            if length(rows) < 1 || length(cols) < 1
                continue
            end
    
            # Move cluster centers to new centroids
            rowmean = clamp(mean(rows), 1, nr)
            colmean = clamp(mean(cols), 1, nc)
            centers[n,:] = round(Int, [rowmean colmean])

            # maxlab[n] = maximum(color_d2[rows,cols]) # slow
        end
        print("Cluster repositioning $iter...")
        toc()
        ########### TODO: put this in a update_centroids!() function ###########
    end
    # newlabels = fixup_clusters(labels, cluster_rows, cluster_cols, npix, centers, s)
    
    # borders = cluster_borders(newlabels)
    # imwrite(grayim(borders), "outliers.jpg")


    # # Visualize original and adjusted grid positions
    # for n = 1:size(grid,1)
    #     i,j = grid[n,:]
    #     img[i,j] = Color.RGB(0,1,0)

    #     # i,j = seeds[n,:]
    #     # img[i,j] = Color.RGB(1,0,0)

    #     i,j = centers[n,:]
    #     img[i,j] = Color.RGB(1,0,1)
    # end

    labels
end

# ###############################################################################
# ###############################################################################
# function fixup_clusters(labels::AbstractArray,
#                         cluster_rows::AbstractArray,
#                         cluster_cols::AbstractArray,
#                         sizes::AbstractArray,
#                         centers::AbstractArray,
#                         min_cluster_size::Integer)
#     nr, nc = size(labels)
#     # newlabels = zeros(Int, nr, nc)
#     newlabels = copy(labels)

#     # Relative neighbor indices
#     i4 = [-1,  0,  1,  0]
#     j4 = [ 0, -1,  0,  1]

#     # newlabels[unlabeled_pixels] = 0

#     # Assign any unassigned pixels to adjacent labels
#     unlabeled_pixels = find(labels .== 0)
#     for p in unlabeled_pixels
#         i,j = ind2sub(size(labels), p)
#         for k = 1:4
#             ni, nj = i+i4[k], j+j4[k]
#             if (1 <= ni <= nr) && (1 <= nj <= nc)
#                 if labels[ni, nj] > 0
#                     newlabels[i,j] = labels[ni, nj]
#                     println("assigning $i,$j to $(labels[ni, nj])")
#                 end
#             end
#         end
#     end
#     unlabeled_pixels = find(labels .== 0)
#     println("$(length(unlabeled_pixels)) unlabeled pixels: $(unlabeled_pixels)")


#     # for n = 1:length(sizes)
#     #     rows = cluster_rows[1:sizes[n], n]
#     #     cols = cluster_cols[1:sizes[n], n]
#     #     # newlabels[rows, cols] = n
#     # end

#     newlabels
# end

# Return row index in the centers array for the point nearest to the nth point.
function nearest_point_index(centers::AbstractArray, n::Integer)
    nn = -1  # Nearest neighbor index
    cn = centers[n,:]
    dmin = convert(Float64, Inf)

    for k = 1:nclusters
        k != n || continue
        r = centers[k,:] - cn
        if dot(r,r) < dmin
            dmin = r
            nn = k
        end
    end
    nn
end

# Assign cluster boundaries at pixels whose label differs from its neighbor 
# below (i+1) or to the right (j+1).
function cluster_borders(labels::AbstractArray)
    m,n = size(labels)
    borders = zeros(m,n)
    for i = 1:m-1
        for j = 1:n-1
            if labels[i,j] != labels[i+1,j] || labels[i,j] != labels[i,j+1] 
                borders[i,j] = 1 
            end
        end
    end
    borders
end

# function regional_adjacency_graph()
#     labels = unique(lpx)
#     nodes = zeros(length(labels), 2)
#     ragmat = copy(boundaries)
#     for label in labels
#         indices = find(lpx .== label)
#         rows, cols = ind2sub(size(lpx), indices)
#         r,c = iround(mean(rows)), iround(mean(cols))
#         nodes[label+1,1:2] = [r c]
#         ragmat[r,c] = 1
#     end
#     imwrite(ragmat, "graph.jpg")
# end

function main()
    img = imread("clutter.jpg")

    # Convert to CIELAB color space for improved gradients and color distances
    imlab = convert(Image{Color.LAB}, float32(img))

    # Approximage number of requested superpixels.
    k = 100

    # Cluster compactness parameter. 
    # Large m favors roundness and uniformity (hex cells as m -> ∞)
    # Small m favors color edge adherence.
    m = 100

    # Note: k and m are coupled, since the clustering distance is a quadrature
    # sum of spatial_distance/s and color_distance/m where s = sqrt(h*w/k).
    
    @time labels = slic(imlab, k, m)

    min_cluster_size = round(Int, sqrt(length(img)/k))
    newlabels, nlabels = reconnect(labels, min_cluster_size)
    
    # show(labels[1:50, 1:50])
    # println("\n nlabels: $nlabels")
    # show(newlabels[1:50, 1:50])

    # Create and save an image
    cluster_img = grayim(newlabels)
    cluster_img /= maximum(cluster_img)
    sc(cluster_img)
    imwrite(cluster_img, "newlabels_grayscale.jpg")

    borders = cluster_borders(labels)
    newborders = cluster_borders(newlabels)

    imwrite(img, "img.jpg")
    imwrite(grayim(borders), "borders.jpg")
    imwrite(grayim(newborders), "newborders.jpg")
end

main()
