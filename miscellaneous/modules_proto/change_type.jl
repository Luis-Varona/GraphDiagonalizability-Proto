function eigvecs_zerooneneg( # Edit docstring
    L::AbstractMatrix{Int64},
    λ_counts::Vector{Pair{Int64, Int64}}
)
    n = size(L, 2)
    eigvecs = [Vector{Int64}[] for _ in λ_counts]
    idx_zero = findfirst(x -> x[1] == 0, λ_counts)
    
    if λ_counts[idx_zero][2] == 1 # The kernel is fully spanned by the all-ones vector
        eigvecs[idx_zero] = [ones(Int64, n)]
    else # The kernel is multi-dimensional, so check for more {-1,0,1}-eigenvectors
        for v in _kernel_eigvecs(n) # All {-1,0,1}-vectors, unique up to span
            iszero(L*v) && push!(eigvecs[idx_zero], v)
        end
    end
    
    # All {-1,0,1}-vectors orthogonal to the all-ones kernel vector, unique up to span
    for v in _laplacian_eigvecs(n)
        # Against all nonzero eigenvalues, check whether `v` is an eigenvector of `L`
        for (idx, λ) in filter(x -> x[2] != 0, collect(enumerate(first.(λ_counts))))
            (L*v == λ*v) && (push!(eigvecs[idx], v); break)
        end
    end
    
    if n % 2 == 0 # Only K_1 and even-ordered graphs have non-kernel {-1,1}-eigenvectors
        idxs_oneneg = findall.(v -> !any(v .== 0), eigvecs)
    else
        idxs_oneneg = [Int64[] for _ in λ_counts]
        idxs_oneneg[idx_zero] = findall(v -> !any(v .== 0), eigvecs[idx_zero])
    end
    
    # Cast each eigenspace's array of eigenvectors to a matrix
    eigvecs = map(vecs -> hcat(vecs...), eigvecs)
    for idx in findall(isempty, eigvecs)
        eigvecs[idx] = Matrix{Int64}(undef, n, 0)
    end
    
    # Check each eigenspace for linearly independent spanning {-1,0,1}- and {-1,1}-sets
    zerooneneg_diag = rank.(eigvecs, 1e-5) .== last.(λ_counts)
    oneneg_diag = rank.(getindex.(eigvecs, :, idxs_oneneg), 1e-5) .== last.(λ_counts)
    
    return (eigvecs, idxs_oneneg, zerooneneg_diag, oneneg_diag)
end