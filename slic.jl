using Images

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

function reconnect(labels::AbstractArray, min_cluster_size::Integer)
    newlabels = zeros(labels)
    nr, nc = size(labels)

    # Relative neighbor indices
    i4 = [-1,  0,  1,  0]
    j4 = [ 0, -1,  0,  1]

    iseg = zeros(Int, nr*nc)
    jseg = zeros(Int, nr*nc)

    label = 1
    adjlabel = 1

    for i = 1:nr
        for j = 1:nc

            # Only enter the loop body if we are starting a new segment
            newlabels[i,j] == 0 || continue
            newlabels[i,j] = label

            # Find an adjacent newlabel for later merging of small segments
            for k = 1:4
                ni, nj = i+i4[k], j+j4[k]
                (1 <= ni <= nr && 1 <= nj <= nc) || continue
                newlabels[ni, nj] > 0 || continue
                adjlabel = newlabels[ni, nj]
            end

            # Do a depth-first search to find connected pixel labels
            iseg[1], jseg[1] = i, j
            counter, segsize = 0, 1
            while counter < segsize
                counter += 1
                for k = 1:4
                    ni, nj = iseg[counter]+i4[k], jseg[counter]+j4[k]

                    (1 <= ni <= nr && 1 <= nj <= nc) || continue
                    labels[i,j] == labels[ni,nj] || continue
                    newlabels[ni,nj] == 0 || continue

                    segsize += 1
                    iseg[segsize], jseg[segsize] = ni, nj
                    newlabels[ni, nj] = label
                end
            end

            # Reassign small clusters to adjlabel
            if segsize < min_cluster_size
                for c = 1:segsize
                    newlabels[iseg[c],jseg[c]] = adjlabel
                end
            else
                label += 1
            end
        end
    end

    newlabels, label
end

function update_distances!(img::AbstractArray, 
                           labels::AbstractArray,
                           centers::AbstractArray,
                           d::AbstractArray,
                           k::Integer,
                           m::Integer)
    nr, nc = size(img)
    nclusters = size(centers, 1)

    # Half-size of clustering window
    s = round(Int, sqrt(nr*nc/k)) - 2
    area = s*s
    m2 = m*m

    # Squared pixel-to-centroid Euclidean distance in LAB color space
    # color_d2 = zeros(Float64, nr, nc)
    # maxlab = 10*10*ones(Float64, nclusters) # Adaptive m-like parameter

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
                end
            end
        end
    end
end

function update_centroids!(labels::AbstractArray, centers::AbstractArray)
    nr, nc = size(labels)
    nclusters = size(centers, 1)
    isums = zeros(Float64, nclusters)
    jsums = zeros(Float64, nclusters)
    npix  = zeros(Int, nclusters)

    # Relative neighbor indices
    i4 = [-1,  0,  1,  0]
    j4 = [ 0, -1,  0,  1]

    for i = 1:nr
        for j = 1:nc
            
            # If pixel was not labeled during the distance update step, try to 
            # assign it an adjacent label.
            if labels[i,j] == 0
                for k = 1:4
                    ni, nj = i+i4[k], j+j4[k]
                    if (1 <= ni <= nr) && (1 <= nj <= nc)
                        if labels[ni, nj] > 0
                            labels[i,j] = labels[ni, nj]
                        end
                    end
                end
            end

            n = labels[i,j]

            n > 0 || continue

            isums[n] += i
            jsums[n] += j
            npix[n] += 1
        end
    end

    for n = 1:nclusters
        npix[n] > 0 || continue
        rowmean = clamp(isums[n]/npix[n], 1, nr)
        colmean = clamp(jsums[n]/npix[n], 1, nc)
        centers[n,:] = round(Int, [rowmean colmean])
    end
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
    
    # Distances from pixels to cluster centers, initialized to ∞
    d = typemax(Float64)*ones(Float64, nr, nc)

    print("Before loop...")
    toc()

    for iter = 1:5

        d *= typemax(Float64)

        tic()
        update_distances!(img, labels, centers, d, k, m)
        print("Dist calcs $iter...")
        toc()

        # TODO: this could go in update_centroids() (?)
        #     # maxlab[n] = maximum(color_d2[rows,cols]) # slow
        tic()
        update_centroids!(labels, centers)
        print("Cluster repositioning $iter...")
        toc()
    end
    
    # Fix any disconnected pixels and merge small clusters to neighbors
    min_cluster_size = round(Int, 0.5*nr*nc/k)
    newlabels, nlabels = reconnect(labels, min_cluster_size)

    newlabels, nlabels
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

# Assign cluster boundaries at pixels whose label differs from its neighbor 
# below (i+1) or to the right (j+1).
function cluster_centroids(labels::AbstractArray, nclusters::Integer)
    nr,nc = size(labels)
    isums = zeros(Float64, nclusters)
    jsums = zeros(Float64, nclusters)
    npix  = zeros(Int, nclusters)
    ctrs  = zeros(Int, nclusters, 2)

    for i = 1:nr
        for j = 1:nc

            n = labels[i,j]

            n > 0 || continue

            isums[n] += i
            jsums[n] += j
            npix[n] += 1
        end
    end

    for n = 1:nclusters
        npix[n] > 0 || continue
        rowmean = clamp(isums[n]/npix[n], 1, nr)
        colmean = clamp(jsums[n]/npix[n], 1, nc)
        ctrs[n,:] = round(Int, [rowmean colmean])
    end
    
    ctrs
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

    # Approximate number of requested superpixels.
    k = 1000

    # Cluster compactness parameter. 
    # Large m favors roundness and uniformity (hex cells as m -> ∞)
    # Small m favors color edge adherence.
    # Note: k and m are coupled, since the clustering distance is a quadrature
    # sum of spatial_distance/s and color_distance/m where s = sqrt(h*w/k).
    m = 10
    
    @time labels, nlabels = slic(imlab, k, m)
    
    # Overlay cluster boundaries on image
    nr, nc = size(labels)
    borders = cluster_borders(labels)
    centroids = cluster_centroids(labels, nlabels)

    # Display cluster boundaries and centroids
    centroid_img = zeros(labels)
    for c in 1:nlabels
        row, col = centroids[c, 1], centroids[c, 2]
        if (1 <= row <= nr) && (1 <= col <= nc)
            centroid_img[row, col] = 1
        end
    end
    segs = Overlay((borders', centroid_img'), 
                   (Color.RGB(0.2,0.2,0.2), Color.RGB(1,1,0)),
                   ((0,1), (0,1))
                   )
                   # (Images.Clamp{Float64}(), Images.Clamp{Float64}()))


    # Superimpose cluster boundaries on image
    for i = 1:nr
        for j = 1:nc
            if borders[i,j] > 0
            img[i,j] = Color.RGB(1,1,1)
            end
        end
    end
    # for c in 1:nlabels
    #     row, col = centroids[c, 1], centroids[c, 2]
    #     if (1 <= row <= nr) && (1 <= col <= nc)
    #         segs[row, col] = Color.RGB(0,1,0)
    #         # centroid_img[row, col] = 1
    #     end
    # end


    imwrite(segs, "segs.jpg")
    imwrite(img, "img.jpg")
    imwrite(grayim(borders), "borders.jpg")
end

main()
