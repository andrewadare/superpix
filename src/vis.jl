function graph_image(g::Graph, centroids::AbstractArray, nr::Int, nc::Int)
    img = zeros(nr, nc)

    # Draw a line for each graph edge
    for e in edges(g)
        r1, c1 = centroids[src(e), 1], centroids[src(e), 2]
        r2, c2 = centroids[dst(e), 1], centroids[dst(e), 2]
        if (1 <= r1 <= nr) && (1 <= c1 <= nc) && (1 <= r2 <= nr) && (1 <= c2 <= nc)
            dist = hypot(r2-r1, c2-c1)
            theta = atan2(r2-r1, c2-c1)
            i, d = 0, 0
            while d < dist
                dr = i*sin(theta)
                dc = i*cos(theta)
                img[round(Int, r1+dr), round(Int, c1+dc)] = 1
                d = hypot(dr, dc)
                i += 1
            end
        end
    end

    # Make a 1-pixel dot at each node
    for c in 1:size(centroids, 1)
        row, col = centroids[c, 1], centroids[c, 2]
        if (1 <= row <= nr) && (1 <= col <= nc)
            img[row, col] = 1
        end
    end

    img
end

function segment_borders(labels::AbstractArray, img::AbstractArray)

    nr, nc = size(img)
    # Relative neighbor indices
    i4 = [-1,  0,  1,  0]
    j4 = [ 0, -1,  0,  1]

    borders = similar(img)
    black = convert(eltype(img), Color.RGB(0,0,0))

    for i = 1:nr
        for j = 1:nc
            borders[i,j] = black
        end
    end

    for i = 1:nr
        for j = 1:nc
            for k = 1:4
                ni, nj = i+i4[k], j+j4[k]

                (1 <= ni <= nr) && (1 <= nj <= nc) || continue

                l, n = labels[i,j], labels[ni,nj]

                if l != n && l > 0
                    borders[i,j] = img[i,j]
                end

            end
        end
    end

    borders
end

function segment_overlay(labels::AbstractArray, 
                         img::AbstractArray,
                         border_colors::AbstractArray,
                         centroids::AbstractArray)

    nr, nc = size(img)
    # Relative neighbor indices
    i4 = [-1,  0,  1,  0]
    j4 = [ 0, -1,  0,  1]

    newimg = copy(img)
    black = convert(eltype(img), Color.RGB(0,0,0))

    # for i = 1:nr
    #     for j = 1:nc
    #         borders[i,j] = black
    #     end
    # end

    for i = 1:nr
        for j = 1:nc
            for k = 1:4
                ni, nj = i+i4[k], j+j4[k]

                (1 <= ni <= nr) && (1 <= nj <= nc) || continue

                l, n = labels[i,j], labels[ni,nj]

                if l != n && l > 0
                    newimg[i,j] = black #border_colors[i,j]
                end

            end
        end
    end

    # Put a crosshair over the segment centroid
    for i = 1:size(centroids, 1)
        r,c = centroids[i,:]
        if r > 0 && c > 0
            labels[r,c] > 0 || continue
            newimg[r,c] = black
            for k = 1:4
                ni, nj = r+i4[k], c+j4[k]
                (1 <= ni <= nr) && (1 <= nj <= nc) || continue
                newimg[ni,nj] = black

                ni, nj = r+2*i4[k], c+2*j4[k]
                (1 <= ni <= nr) && (1 <= nj <= nc) || continue
                newimg[ni,nj] = black

            end

        end
    end

    newimg
end
