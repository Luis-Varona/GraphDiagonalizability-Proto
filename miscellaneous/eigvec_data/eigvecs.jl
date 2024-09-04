include("../../modules/LaplacianSpectra.jl")
laplacian_eigvecs = LaplacianSpectra.laplacian_eigvecs

using CodecZlib
using JLD2

function main()
    for n in 1:20
        E = hcat(Vector{Vector{Int64}}(collect(laplacian_eigvecs(n)))...)
        jldsave("miscellaneous/eigvec_data/eigvecs_$n.jld2", true; E)
        println("Finished n = $n")
    end
end