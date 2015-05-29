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

            # Depth-first search: find all connected labels matching labels[i,j].
            # The pixel rows and columns of this cluster/segment will be stored
            # in iseg[1:segsize] and jseg[1:segsize].
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

    newlabels, nlabels = slic(img, centers, k, m, niter=5)
end

function slic(img::AbstractArray, 
              seeds::AbstractArray,
              k::Integer,
              m::Integer;
              niter = 5)
    nr, nc    = size(img)
    nclusters = size(seeds, 1)

    # Cluster labels for each pixel
    labels = zeros(Int, nr, nc)
    
    # Distances from pixels to cluster seeds, initialized to ∞
    d = typemax(Float64)*ones(Float64, nr, nc)

    for iter = 1:niter

        d *= typemax(Float64)

        update_distances!(img, labels, seeds, d, k, m)
        update_centroids!(labels, seeds)
    end
    
    # Fix any disconnected pixels and merge small clusters to neighbors
    min_cluster_size = round(Int, 0.5*nr*nc/k)
    newlabels, nlabels = reconnect(labels, min_cluster_size)

    newlabels, nlabels
end
