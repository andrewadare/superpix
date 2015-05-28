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

function lab_array(rgb_img_name::AbstractString, rows::UnitRange, cols::UnitRange)
    img = subim(imread(rgb_img_name), "x", cols, "y", rows)
    nr, nc = size(img)
    alab = Array(Color.Lab{Float32}, nr, nc)

    for i in nr
        for j in nc
            alab[i,j] = convert(Color.Lab{Float32}, img[i,j])
        end
    end

    alab'
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
