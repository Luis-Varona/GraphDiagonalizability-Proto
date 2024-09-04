using BenchmarkTools

include("../modules_master/LaplacianSpectra.jl")
using .LaplacianSpectra: eigvecs_zerooneneg, integer_eigvals_by_multi

# function main()
#     for l in first.(M[3])
#         eigvecs_zerooneneg(L, l)
#     end
# end

L = [[ 5 -1 -1 -1 -1 -1  0  0  0  0  0  0]
[-1  5 -1  0  0  0 -1  0 -1 -1  0  0]
[-1 -1  5  0  0  0  0 -1  0  0 -1 -1]
[-1  0  0  5 -1 -1 -1 -1  0  0  0  0]
[-1  0  0 -1  5 -1  0  0 -1  0 -1  0]
[-1  0  0 -1 -1  5  0  0  0 -1  0 -1]
[ 0 -1  0 -1  0  0  5 -1 -1 -1  0  0]
[ 0  0 -1 -1  0  0 -1  5  0  0 -1 -1]
[ 0 -1  0  0 -1  0 -1  0  5 -1 -1  0]
[ 0 -1  0  0  0 -1 -1  0 -1  5  0 -1]
[ 0  0 -1  0 -1  0  0 -1 -1  0  5 -1]
[ 0  0 -1  0  0 -1  0 -1  0 -1 -1  5]]

M = integer_eigvals_by_multi(L)[3]

# function main()
#     eigvecs = Vector{Vector{Vector{Int64}}}(undef, length(M))
#     eigvals = Vector{Vector{Int64}}(undef, length(M))
    
#     for (i, l) in enumerate(first.(M))
#         (eigvecs[i], eigvals[i]) = eigvecs_zerooneneg(L, l)
#     end
    
#     return (eigvecs, eigvals)
# end

p = [0 0 0 0; 0 0 0 0; 0 0 1 -1; 0 0 -1 1]
P = integer_eigvals_by_multi(p)[3]