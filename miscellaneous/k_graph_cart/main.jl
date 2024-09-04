#- Imports
include("../../modules/KOrthogonalizability.jl")
is_quasi_orthogonalizable = KOrthogonalizability.is_quasi_orthogonalizable

using CodecZlib
using Combinatorics: combinations, permutations, multiset_permutations
using JLD2: jldsave, load
using LinearAlgebra: I, rank


#- Main function
function main()::Nothing
    n = 5
    L = complete_cart(n)
    μ = (n - 1)^2
    
    # E = get_2n_eigvecs(n)
    # save_2n_eigvecs(n, "miscellaneous/k_graph_cart/eigvecs.jld2") ###
    E = load("miscellaneous/k_graph_cart/eigvecs.jld2")["E"] ###
    E = Vector{Int64}.(eachcol(E)) ###
    
    (is_WHD, basis) = is_WHD_eigenspace(L, μ, E)
    println(is_WHD)
    jldsave("basis.jld2"; basis)
end


#- Helper functions
kron(a...) = Base.kron(a...)
function kron_sum(A::Matrix{Int64}, B::Matrix{Int64})::Matrix{Int64}
    return kron(A, I(size(B, 1))) + kron(I(size(A, 1)), B)
end

k_graph(n::Int64)::Matrix{Int64} = n*I(n) - ones(Int64, n, n)
complete_cart(n::Int64)::Matrix{Int64} = kron_sum(k_graph(n), k_graph(n))


#-
function all_laplacian_eigvecs(n::Int64)::Iterators.Flatten
    function entries_set(k::Int64)
        entries = Vector{Int64}(undef, n)
        entries[1:k] .= 1
        entries[(k + 1):2k] .= -1
        entries[(2k + 1):n] .= 0
        return entries
    end
    
    return Iterators.flatten(
        multiset_permutations.(entries_set.(0:floor(Int64, n/2)), n)
    )
end


#-
function get_2n_eigvecs(n::Int64)::Vector{Vector{Int64}}
    eigvec_blocks = collect(all_laplacian_eigvecs(n))
    eigvecs = Vector{Int64}[]
    
    for mat in Iterators.product(Iterators.repeated(eigvec_blocks, n)...)
        if iszero(sum(reduce(hcat, mat), dims = 2))
            v = vec(reduce(hcat, mat))
            if !iszero(v) && v[findfirst(x -> x != 0, v)] == 1
                push!(eigvecs, v)
            end
        end
    end
    
    return eigvecs
end


#-
function is_WHD_eigenspace(
    L::AbstractMatrix{Int64},
    μ::Int64,
    eigvecs::Vector{Vector{Int64}}
)::Tuple{Bool, Matrix{Int64}}
    
    num_eigvecs = length(eigvecs)
    
    function DFS(idxs::Vector{Int64})::Tuple{Bool, Vector{Int64}}
        depth = length(idxs)
        eigbasis = hcat(eigvecs[idxs]...) # The current eigenbasis submatrix
        
        if depth == 10 ###
            println(idxs) ###
        end ###
        
        if rank(eigbasis, 1e-5) != depth # Verify linear independence
            is_quasi = false
            idxs_final = Int64[]
        else
            (is_quasi, order) = is_quasi_orthogonalizable(eigbasis)
            
            if !is_quasi
                idxs_final = Int64[]
            elseif depth == μ # A complete set of basis vectors has been found
                idxs_final = idxs[order] # Impose a k-orthogonal column ordering
            else
                is_quasi = false
                idxs_final = Int64[]
                last = idxs[end] # The greatest index in the (sorted) list
                
                # Only add subsequent indices to avoid duplicate combinations
                for idx_next in (last + 1):num_eigvecs
                    idxs_next = vcat(idxs, idx_next)
                    
                    # Recursively call `DFS` and search for a complete basis set
                    (is_quasi, idxs_final) = DFS(idxs_next)
                    if is_quasi
                        break
                    end
                end
            end
        end
        
        return (is_quasi, idxs_final)
    end
    
    is_WHD = false
    
    for root in combinations(1:num_eigvecs, 2)
        (is_WHD, idxs_final) = DFS(root)
        if is_WHD
            basis = hcat(eigvecs[idxs_final]...)
            break
        end
    end
    
    if !is_WHD
        basis = Matrix{Int64}(undef, size(L, 1), 0)
    end
    
    return (is_WHD, basis)
end


#-
function save_2n_eigvecs(n::Int64, path::String)::Nothing
    isdir(path) || mkpath(path)
    E = get_2n_eigvecs(n)
    E = hcat(E...)
    jldsave(path, true; E)
end


#- Run the main function
# main()