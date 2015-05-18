using Images

include("/Users/adare/repos/superpix/slic.jl")

function main()
    a = readdlm("pbuff.txt", UInt32)
    nrows, ncols = 171, 361
    a = reshape(a, ncols, nrows)
    img = reinterpret(Color.RGB24, Image(a))
    

    show(img)
    println("")    
    imsep = separate(img)
    show(imsep)
    # imlab = convert(Image{Color.LAB}, map(Float32, separate(img))) # No method
    # if spatialorder(imlab) != spatialorder(img)
    #     perm = spatialpermutation(spatialorder(img), imlab)
    #     imlab = permutedims(imlab, perm)
    # end
    # @assert spatialorder(img) == spatialorder(imlab)
    return










    alab = Array(Color.Lab{Float32}, ncols, nrows)
    
    for i = 1:ncols
        for j = 1:nrows
            alab[i,j] = convert(Color.Lab{Float32}, img[i,j])
        end
    end
    
    imlab = Image(alab)
    imlab["IMcs"] = "sRGB"
    imlab["spatialorder"] = ["x","y"]
    imlab["pixelspacing"] = [1,1]
    show(imlab)

    k, m = 1000, 10
    @time labels, nlabels = slic(imlab, k, m)
    
    # imwrite(imlab, "imlab.jpg")
end

main()
