using Images

depth_img = FileIO.load("kinect-depth.png")
color_img = FileIO.load("kinect-color.png")

# Relative neighbor indices
i4 = [-1,  0,  1,  0]
j4 = [ 0, -1,  0,  1]

# Convert depth image to a plain old array to simplify access
dep_x_range = 60:420
dep_y_range = 120:290
dep = convert(Array{Float32}, depth_img[dep_x_range, dep_y_range])
rgb_x_range = 390:1445
rgb_y_range = 250:750
rgb = color_img[rgb_x_range, rgb_y_range]

# Assign empty pixels to a neighboring value
empties = find(dep .< 0.001)
for e in empties
    println(e, " ", dep[e])

    i,j = ind2sub(size(dep), e)

    for k = 1:4
        ni, nj = i+i4[k], j+j4[k]
        if (1 <= ni <= size(dep, 1)) && (1 <= nj <= size(dep, 2))
            if dep[ni, nj] > 0
                dep[i,j] = dep[ni, nj]
            end
        end
    end

end

# Adjust grayscale range to 0,1
dep -= minimum(dep)
dep /= maximum(dep)

println("cropped depth image aspect ratio: $(size(dep, 1)/size(dep, 2))")
println("cropped rgb   image aspect ratio: $(size(rgb, 1)/size(rgb, 2))")

imwrite(dep', "dep.png") # Size is 171x361
imwrite(rgb', "rgb.png") # Size is 501x1056 -> Resize to 171x361
