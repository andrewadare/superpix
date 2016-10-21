# Read in color image and convert to CIE Lab color space
function lab_image(imgname::AbstractString, rows::UnitRange, cols::UnitRange)
    imrgb = FileIO.load(imgname)
    imlab = convert(Image{Colors.Lab}, imrgb["x", cols, "y", rows])
    # imlab = convert(Image{Colors.Lab}, map(Float32, separate(imrgb["x", cols, "y", rows])))
    imlab.properties["spatialorder"] = ["x","y"]
    imlab
end

# Read in depth image and preprocess to improve contrast
function depth_image(imgname::AbstractString, rows::UnitRange, cols::UnitRange)
    imdep = FileIO.load(imgname)
    dep = convert(Array, map(Float32, separate(imdep["x", cols, "y", rows])))

    # Remove singleton dimension
    dep = reshape(dep, length(cols), length(rows))

    # println("Before rescaling: min, max = ", minimum(dep), " ", maximum(dep))

    for i=1:2
        mindep = minimum(dep[dep .> 0])
        dep = clamp(dep, mindep, maximum(dep)) - mindep
        dep /= maximum(dep)
    end

    # Gaussian blur helps to smooth over empty pixels and noise
    dep = imfilter_gaussian(dep, [3,3])
    mindep = minimum(dep[dep .> 0])
    dep = clamp(dep, mindep, maximum(dep)) - mindep
    dep /= maximum(dep)
    # println("After rescaling:  min, max = ", minimum(dep), " ", maximum(dep))
    dep
end

function unpack_drgb_array(pbuff::Array{UInt32, 1}, nrows::Int64, ncols::Int64)

    pbuff_2d = reshape(pbuff, ncols, nrows)

    # Unpack UInt32s into LAB color and depth arrays.
    dep  = Array(Float32, ncols, nrows)
    alab = Array(Colors.Lab{Float32}, ncols, nrows)

    for i = 1:ncols
        for j = 1:nrows
            p = UInt32(pbuff_2d[i,j])
            b = ((p >>  0) & 0xFF)/255
            g = ((p >>  8) & 0xFF)/255
            r = ((p >> 16) & 0xFF)/255
            d = ((p >> 24) & 0xFF)/255
            dep[i,j] = d
            alab[i,j] = convert(Colors.Lab{Float32}, Colors.RGB(r, g, b))
        end
    end
    dep, alab
end

function lab_array(rgb_img_name::AbstractString, rows::UnitRange, cols::UnitRange)
    img = subim(load(rgb_img_name), "x", cols, "y", rows)
    nr, nc = size(img)
    alab = Array(Colors.Lab{Float32}, nr, nc)

    for i in nr
        for j in nc
            alab[i,j] = convert(Colors.Lab{Float32}, img[i,j])
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
