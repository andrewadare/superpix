
# Compute depth image gradients
function depth_image_grad(dep::AbstractArray)
    gx, gy = imgradients(dep, "ando5")
    depgrad = magnitude(data(gx),data(gy))

    # Scale gradient magnitudes to range from 0 to 1
    mingrad = minimum(depgrad[depgrad .> 0])
    depgrad = clamp(depgrad, mingrad, maximum(depgrad)) - mingrad
    depgrad /= maximum(depgrad)
    depgrad
end

# Compute color image gradient (mean over L,a,b magnitudes)
function color_image_grad(imlab::AbstractArray)
    gx, gy = imgradients(imlab, "ando5")
    gradmag = magnitude(data(gx),data(gy))
    gradmean = reshape(mean(gradmag, colordim(gx)), tuple(size_spatial(gx)...))

    # Scale gradient magnitudes to range from 0 to 1
    mingrad = minimum(gradmean[gradmean .> 0])
    gradmean = clamp(gradmean, mingrad, maximum(gradmean)) - mingrad
    gradmean /= maximum(gradmean)
    gradmean
end
