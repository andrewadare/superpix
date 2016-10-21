function generalized_mean{T<:Images.ColorTypes.Lab, N}(A::AbstractArray{T, N})
    c = zeros(Float32, 3)
    for pixel in A
        c[1] += pixel.l
        c[2] += pixel.a
        c[3] += pixel.b
    end
    c /= length(A)
    return Colors.Lab{Float32}(c[1], c[2], c[3])
end

# function generalized_mean{T<:Union(Number, Images.ColorTypes.Gray, Images.ColorTypes.Rgb)}(A::AbstractArray{T, N})
function generalized_mean(A::AbstractArray)
    # t = eltype(A)
    # if t <: Union(Number, Images.ColorTypes.Gray, Images.ColorTypes.RGB)
        return mean(A)
    # elseif t <: Images.ColorTypes.Lab
        # c = zeros(Float32, 3)
        # for pixel in A
        #     c[1] += pixel.l
        #     c[2] += pixel.a
        #     c[3] += pixel.b
        # end
        # c /= length(A)
        # return Colors.Lab{Float32}(c[1], c[2], c[3])
    # end

    # error("No implementation for type $t")
    # NaN
end

function color_means(img, labels, nlabels)
    # d = Dict{Int, Float32}()
    d = Dict()
    superpx_mean_img = similar(img)

    var = zeros(Float64, 3)
    for label in 1:nlabels
        indices = find(labels .== label)
        segment = img[indices]
        mu = generalized_mean(segment)
        # println(mu)

        d[label] = mu
        superpx_mean_img[indices] = mu
    end
    d, superpx_mean_img
end

function color_moments(img, labels, nlabels)
    superpx_mean_img = zeros(img)
    superpx_std_img = zeros(img)

    var = zeros(Float64, 3)
    for label in 1:nlabels
        indices = find(labels .== label)
        segment = img[indices]
        mu = mean(segment)
        # sd = std(segment) # No implementation for ColorValues

        isnan(mu) && continue

        # Compute the variance in color over this superpixel.
        # The color variance is a 3-component vector of floats for R,G,B.
        cdiff = Array(Float64, length(segment), 3)
        for p = 1:length(segment)
            c = img[p]
            cdiff[p, :] = Float64[c.r - mu.r c.g - mu.g c.b - mu.b]
        end
        for color = 1:3
            var[color] = mapreduce(x->x^2, +, cdiff[:, color])/max(1, length(segment)-1)
        end

        stdev = sqrt(var)
        superpx_mean_img[indices] = mu
        superpx_std_img[indices] = Colors.RGB(stdev[1], stdev[2], stdev[3])
    end
    superpx_mean_img, superpx_std_img
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

