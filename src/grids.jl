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

# Return grid points nudged to the minimum-gradient position in a 3x3 window
function adjusted_grid(grid::AbstractArray, grad_mag::AbstractArray)

    nr, nc = size(grad_mag)
    newgrid = copy(grid)

    for n = 1:size(grid,1)

        # Make a 3x3 window around the grid point i,j.
        i,j = grid[n,:]
        rows = max(1,i-1) : min(nr,i+1)
        cols = max(1,j-1) : min(nc,j+1)
        window = view(grad_mag, rows, cols)

        # Find pixel position in window with the smallest gradient.
        wi, wj = ind2sub(size(window), indmin(window))

        newgrid[n,:] = [rows[1] - 1 + wi, cols[1] - 1 + wj]
    end

    newgrid
end

# Compute squared magnitude of image gradient for use in adjusted_grid
function gradmag2(img::AbstractArray, method::String)
    gx, gy = imgradients(img, method)
    g = magnitude(data(gx), data(gy))

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
