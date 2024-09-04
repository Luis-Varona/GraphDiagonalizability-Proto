#-
include("../../modules/KOrthogonalizability.jl")
using .KOrthogonalizability: is_k_orthogonalizable

using Combinatorics
using JLD2
using Kronecker
using LinearAlgebra

#-
E = load("miscellaneous/k_graph_cart/eigvecs.jld2")["E"]
E = Vector{Int64}.(eachcol(E))

#-
function stdBasisVec(i::Int64, n::Int64)
    e = zeros(Int64, n)
    e[i] = 1
    return e
end

#-
function NJ_vec(i::Int64, j::Int64, n::Int64)
    arg1 = stdBasisVec(i, n) - stdBasisVec(i + 1, n)
    arg2 = stdBasisVec(j, n) - stdBasisVec(j + 1, n)
    return arg1 ⊗ arg2
end

#-
function get_basis(n::Int64)
    V = Matrix{Int64}(undef, n^2, (n - 1)^2)
    
    for i in 1:(n - 1)
        for j in 1:(n - 1)
            V[:, (i - 1)*(n - 1) + j] = NJ_vec(i, j, n)
        end
    end
    
    return V
end

#-

#-
# for n = 1:7
#     println("n = ", n, ": ", is_quasi_orthogonalizable(get_basis(n)))
# end