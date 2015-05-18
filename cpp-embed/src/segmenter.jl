using Images

include("/Users/adare/repos/superpix/slic.jl")

function main(a::Array{UInt32, 1}, nrows::Int64, ncols::Int64)
    a = reshape(a, ncols, nrows)'
    img = reinterpret(Color.RGB24, Image(a))

    # Convert to CIELAB color space for improved gradients and color distances.
    # Correct back to x vs y order of img if not already matching.
    imlab = convert(Image{Color.LAB}, map(Float32, separate(img)))
    if spatialorder(imlab) != spatialorder(img)
        perm = spatialpermutation(spatialorder(img), imlab)
        imlab = permutedims(imlab, perm)
    end
    @assert spatialorder(img) == spatialorder(imlab)

    imwrite(imlab, "imlab.jpg")
end
