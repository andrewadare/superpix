using Images

include("/Users/adare/repos/superpix/slic.jl")

function main()
    a = readdlm("../pbuff.txt", UInt32)
    nrows, ncols = 171, 361
    a = reshape(a, ncols, nrows)
    # img = reinterpret(Color.RGB24, Image(a))
    # img = reinterpret(Color.ARGB32, Image(a))
    
    # show(img)
    # println("")    
    # imsep = separate(img)
    # show(imsep)
    # imlab = convert(Image{Color.LAB}, map(Float32, separate(img))) # No method
    # if spatialorder(imlab) != spatialorder(img)
    #     perm = spatialpermutation(spatialorder(img), imlab)
    #     imlab = permutedims(imlab, perm)
    # end
    # @assert spatialorder(img) == spatialorder(imlab)
    # return

    alab = Array(Color.Lab{Float32}, ncols, nrows)
    dep  = Array(Float32, ncols, nrows)
    rgb  = Array(Color.RGB, ncols, nrows)

    for i = 1:ncols
        for j = 1:nrows
            p = UInt32(a[i,j])
            b = ((p >>  0) & 0xFF)/255
            g = ((p >>  8) & 0xFF)/255
            r = ((p >> 16) & 0xFF)/255
            d = ((p >> 24) & 0xFF)/255

            # alab[i,j] = convert(Color.Lab{Float32}, Color.RGB(p.r, p.g, p.b))
            rgb[i,j] = Color.RGB(r,g,b)
            dep[i,j] = d
        end
    end
    depth_img = grayim(dep)
    imwrite(depth_img, "depth.jpg")    
    # imlab = Image(alab)
    # imlab["IMcs"] = "sRGB"
    # imlab["spatialorder"] = ["x","y"]
    # imlab["pixelspacing"] = [1,1]
    # show(imlab)

    # k, m = 1000, 10
    # @time labels, nlabels = slic(imlab, k, m)
    
    # imwrite(imlab, "imlab.jpg")
end

main()
