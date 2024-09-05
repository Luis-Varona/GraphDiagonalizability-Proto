#- FUNCTION: `_generalgraph_eigvecs_zerooneneg`
"""
    _generalgraph_eigvecs_zerooneneg(L, λ_counts)

ADD LATER

# Arguments
- `L::AbstractMatrix{Int64}`: ADD LATER
- `λ_counts::Vector{Pair{Int64, Int64}}`: ADD LATER

# Returns
- `eigvecs::Vector{Matrix{Int64}}`: ADD LATER
- `cols_oneneg::Vector{Vector{Int64}}`: ADD LATER
- `zerooneneg_bases::BitVector`: ADD LATER
- `oneneg_bases::BitVector`: ADD LATER
"""
function _generalgraph_eigvecs_zerooneneg(
    L::AbstractMatrix{Int64},
    λ_counts::Vector{Pair{Int64, Int64}},
)
    λs_unique = first.(λ_counts)
    multis = last.(λ_counts)
    
    n = size(L, 1)
    idx_kernel = findfirst(iszero, λs_unique)
    # Use an array of resizeable matrices to store an unknown number of eigenvectors
    eigvecs_elastic = [ElasticMatrix{Int64}(undef, n, 0) for _ in λs_unique]
    
    if multis[idx_kernel] == 1 # The all-ones vector fully spans the kernel
        append!(eigvecs_elastic[idx_kernel], ones(Int64, n))
    else # Find all {-1,0,1}-vectors unique up to span in the multi-dimensional kernel
        for v in potential_kernel_eigvecs(n)
            iszero(L*v) && append!(eigvecs_elastic[idx_kernel], v)
        end
    end
    
    # Check all {-1,0,1}-vectors orthogonal to the all-ones vector, unique up to span
    for v in potential_nonkernel_eigvecs(n)
        # Against all nonzero eigenvalues, check whether `v` is an eigenvector of `L`
        for (idx, λ) in deleteat!(collect(enumerate(λs_unique)), idx_kernel)
            (L*v == λ*v) && (append!(eigvecs_elastic[idx], v); break)
        end
    end
    
    eigvecs = Matrix{Int64}.(eigvecs_elastic) # Convert back to the standard matrix type
    cols_oneneg = _get_cols_oneneg(eigvecs, idx_kernel, n)
    
    # Test for linearly independent spanning {-1,0,1}- and {-1,1}-eigenbases
    zerooneneg_bases = (rank.(eigvecs, 1e-5) .== multis)
    oneneg_bases = (rank.(getindex.(eigvecs, :, cols_oneneg), 1e-5) .== multis)
    
    return (eigvecs, cols_oneneg, zerooneneg_bases, oneneg_bases)
end


#- FUNCTION: `_completegraph_eigvecs_zerooneneg`
"""
    _completegraph_eigvecs_zerooneneg(λ_counts)

ADD LATER

# Arguments
- `λ_counts::Vector{Pair{Int64, Int64}}`: ADD LATER

# Returns
- `eigvecs::Vector{Matrix{Int64}}`: ADD LATER
- `cols_oneneg::Vector{Vector{Int64}}`: ADD LATER
- `zerooneneg_bases::BitVector`: ADD LATER
- `oneneg_bases::BitVector`: ADD LATER
"""
function _completegraph_eigvecs_zerooneneg(λ_counts::Vector{Pair{Int64, Int64}})
    n = sum(last.(λ_counts)) # Every order n graph has exactly n Laplacian eigenvalues
    idx_kernel = findfirst(x -> iszero(x[1]), λ_counts)
    # Complete graph Laplacians have exactly one non-kernel eigenspace for all n ≥ 2
    idx_nonkernel = (idx_kernel == 1) ? 2 : 1
    
    eigvecs = Vector{Matrix{Int64}}(undef, 2)
    eigvecs[idx_kernel] = ones(Int64, n, 1) # The all-ones vector fully spans the kernel
    # All vectors orthogonal to the all-ones kernel vector are non-kernel eigenvectors
    eigvecs[idx_nonkernel] = hcat(collect(potential_nonkernel_eigvecs(n))...)
    
    cols_oneneg = _get_cols_oneneg(eigvecs, idx_kernel, n)
    zerooneneg_bases = trues(2) # Both eigenspaces have {-1,0,1}-bases for all n ≥ 2
    
    if n % 2 == 0 # Both eigenspaces have {-1,1}-bases for all even n ≥ 2
        oneneg_bases = trues(2)
    else # No odd-dimensional {-1,1}-vectors are orthogonal to the all-ones vector
        oneneg_bases = BitVector(undef, 2)
        oneneg_bases[idx_kernel] = true
        oneneg_bases[idx_nonkernel] = false
    end
    
    return (eigvecs, cols_oneneg, zerooneneg_bases, oneneg_bases)
end