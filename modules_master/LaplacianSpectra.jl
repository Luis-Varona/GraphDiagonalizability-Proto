module LaplacianSpectra
    #- EXPORTS AND IMPORTS
    export is_spectrum_integral
    export eigvecs_zerooneneg
    
    include("EigenvectorGenerators.jl")
    using .EigenvectorGenerators: potential_kernel_eigvecs, potential_nonkernel_eigvecs
    
    using Combinatorics: combinations
    using DataStructures: counter
    using ElasticArrays: ElasticMatrix
    using LinearAlgebra: eigvals, rank
    using RowEchelon: rref_with_pivots
    
    
    #- FUNCTION: `is_spectrum_integral`
    """
        is_spectrum_integral(X)
    
    Determine whether an integer matrix is spectrum integral and find its eigenvalues if so.
    
    If an undirected Laplacian matrix has any non-integer eigenvalues, the graph it
    represents cannot be `{-1,0,1}`-diagonalizable (assuming integer edge weights).
    
    # Arguments
    - `X::AbstractMatrix{Int64}`: an entrywise-integer matrix. In the context of the larger
        `GraphDiagonalizability` package, this will typically be the Laplacian matrix of
        some undirected graph with integer edge weights.
    
    # Returns
    - `spectrum_integral::Bool`: whether the eigenvalues of `X` are all integers.
    - `λs_sorted::Vector{Int64}`: if `spectrum_integral`, an array of all eigenvalues of `X`
        sorted by multiplicity then by value. Otherwise, an empty array.
    - `λ_counts::Vector{Pair{Int64, Int64}}`: if `spectrum_integral`, an array of all unique
        eigenvalues of `X` and their multiplicities in the same order as `λs_sorted`.
        Otherwise, an empty array.
    
    # Examples
    Test the unweighted Laplacian matrix of the cycle graph `C_6`:
    ```jldoctest
    julia> L = [ 2  -1  -1   0   0   0;
                -1   2   0   0  -1   0;
                -1   0   2   0   0  -1;
                 0   0   0   2  -1  -1;
                 0  -1   0  -1   2   0;
                 0   0  -1  -1   0   2];
    
    julia> is_spectrum_integral(L)
    (true, [0, 4, 1, 1, 3, 3], [0 => 1, 4 => 1, 1 => 2, 3 => 2])
    ```
    
    Test a weighted Laplacian matrix of the star graph `S_5`:
    ```jldoctest
    julia> L = [ 1   0   0   0  -1;
                 0   1   0   0  -1;
                 0   0   1   0  -1;
                 0   0   0   7  -7;
                -1  -1  -1  -7  10];
    
    julia> is_spectrum_integral(L)
    (false, Int64[], Pair{Int64, Int64}[])
    ```
    
    Test a weighted Laplacian matrix of the strong product `K_2 ⊠ K_{1,1,2}`:
    ```jldoctest
    julia> L = [ 6  -2   0   0  -1  -1  -1  -1;
                -2   6   0   0  -1  -1  -1  -1;
                 0   0   6  -2  -1  -1  -1  -1;
                 0   0  -2   6  -1  -1  -1  -1;
                -1  -1  -1  -1   9  -3  -1  -1;
                -1  -1  -1  -1  -3   9  -1  -1;
                -1  -1  -1  -1  -1  -1   9  -3;
                -1  -1  -1  -1  -1  -1  -3   9];
    
    julia> is_spectrum_integral(L)
    (true, [0, 4, 12, 12, 8, 8, 8, 8], [0 => 1, 4 => 1, 12 => 2, 8 => 4])
    ```
    """
    function is_spectrum_integral(X::AbstractMatrix{Int64})
        # Eigenvalue computations are imprecise, so check whether they are close to integers
        λs_float = eigvals(X)
        λs_int = round.(Int64, real.(λs_float))
        spectrum_integral = isapprox(λs_float, λs_int)
        
        if !spectrum_integral
            λs_sorted = Int64[]
            λ_counts = Pair{Int64, Int64}[]
        else
            # Sort all unique eigenvalues by ascending multiplicity, then by ascending value
            λ_counts = sort(sort(collect(counter(λs_int))), by = last)
            # Sort all eigenvalues, including repetitions, in the same order as `λ_counts`
            λs_sorted = Int64.(vcat(fill.(first.(λ_counts), last.(λ_counts))...))
        end
        
        return (spectrum_integral, λs_sorted, λ_counts)
    end
    
    
    #- FUNCTION: `eigvecs_zerooneneg`
    """
        eigvecs_zerooneneg(L, λ_counts)
    
    Find all normalized `{-1,0,1}`- and `{-1,1}`-eigenvectors of a Laplacian matrix.
    
    The first nonzero entry of each eigenvector is normalized to `1`, enforcing pairwise
    linear independence while still covering the equivalence class of span.
    
    Weighted Laplacian matrices are supported, but all edge weights must be integers.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of some undirected graph.
    - `λ_counts::Vector{Pair{Int64, Int64}}`: the unique eigenvalues of `L` and their
        corresponding multiplicities, sorted by frequency then by value.
    
    # Returns
    - `eigvecs::Vector{Matrix{Int64}}`: an array of all `{-1,0,1}`-eigenvectors (as columns
        of matrices) corresponding to each unique eigenvalue in `λ_counts`.
    - `cols_oneneg::Vector{Vector{Int64}}`: an array of all column indices of each matrix
        in `eigvecs` corresponding to `{-1,1}`-eigenvectors in that eigenspace.
    - `bases_zerooneneg::BitVector`: a Boolean array indicating whether each eigenspace has
        enough linearly independent `{-1,0,1}`-vectors to form a basis set.
    - `bases_oneneg::BitVector`: a Boolean array indicating whether each eigenspace has
        enough linearly independent `{-1,1}`-vectors to form a basis set.
    
    # Examples
    Test the unweighted Laplacian matrix of the cograph `K_1 ⊔ K_1 ⊔ K_3`:
    ```jldoctest
    julia> L = [0  0   0   0   0;
                0  0   0   0   0;
                0  0   2  -1  -1;
                0  0  -1   2  -1
                0  0  -1  -1   2];
    
    julia> λ_counts = is_spectrum_integral(L)[3]
    2-element Vector{Pair{Int64, Int64}}:
     3 => 2
     0 => 3
    
    julia> (eigvecs, cols_oneneg, bases_zerooneneg, bases_oneneg) = eigvecs_zerooneneg(
               L,
               λ_counts,
           );
    
    julia> eigvecs[1]
    5×3 Matrix{Int64}:
      0   0   0
      0   0   0
      1   1   0
     -1   0   1
      0  -1  -1
    
    julia> eigvecs[2]
    5×13 Matrix{Int64}:
     1   1  1   1   1   1  1   1  1  0   0  0  0
     1  -1  0   1  -1   0  1  -1  0  1   1  1  0
     1   1  1  -1  -1  -1  0   0  0  1  -1  0  1
     1   1  1  -1  -1  -1  0   0  0  1  -1  0  1
     1   1  1  -1  -1  -1  0   0  0  1  -1  0  1
    
    julia> (cols_oneneg, bases_zerooneneg, bases_oneneg)
    ([Int64[], [1, 2, 4, 5]], Bool[1, 1], Bool[0, 1])
    ```
    
    Test a weighted Laplacian matrix of the strong product `P_2 ⊠ P_3`:
    ```jldoctest
    julia> L = [ 3  -1   0   0  -1  -1;
                -1   3   0   0  -1  -1;
                 0   0   3  -1  -1  -1;
                 0   0  -1   3  -1  -1;
                -1  -1  -1  -1   6  -2;
                -1  -1  -1  -1  -2   6];
    
    julia> λ_counts = is_spectrum_integral(L)[3]
    5-element Vector{Pair{Int64, Int64}}:
     0 => 1
     2 => 1
     6 => 1
     8 => 1
     4 => 2
    
    julia> (eigvecs, cols_oneneg, bases_zerooneneg, bases_oneneg) = eigvecs_zerooneneg(
               L,
               λ_counts,
           );
    
    julia> eigvecs[1]
    6×1 Matrix{Int64}:
     1
     1
     1
     1
     1
     1
    
    julia> eigvecs[2]
    6×1 Matrix{Int64}:
      1
      1
     -1
     -1
      0
      0
    
    julia> eigvecs[3]
    6×0 Matrix{Int64}
    
    julia> eigvecs[4]
    6×1 Matrix{Int64}:
      0
      0
      0
      0
      1
     -1
    
    julia> eigvecs[5]
    6×4 Matrix{Int64}:
      1   1   1   0
     -1  -1  -1   0
      0  -1   1   1
      0   1  -1  -1
      0   0   0   0
      0   0   0   0
    
    julia> (cols_oneneg, bases_zerooneneg, bases_oneneg)
    ([[1], Int64[], Int64[], Int64[], Int64[]], Bool[1, 1, 0, 1, 1], Bool[1, 0, 0, 0, 0])
    ```
    """
    function eigvecs_zerooneneg(
        L::AbstractMatrix{Int64},
        λ_counts::Vector{Pair{Int64, Int64}},
    )
        n = size(L, 1)
        
        if n == 0 # The Laplacian matrix of the null graph has no eigenspaces
            eigvecs = Matrix{Int64}[]
            cols_oneneg = Vector{Int64}[]
            bases_zerooneneg = trues(0)
            bases_oneneg = trues(0)
        # Every vector is an eigenvector of the zero matrix, which only has a 0-eigenspace
        elseif iszero(L)
            eigvecs = [hcat(collect(potential_kernel_eigvecs(n))...)]
            cols_oneneg = [findall(v -> !any(v .== 0), eachcol(eigvecs[1]))]
            bases_zerooneneg = trues(1)
            bases_oneneg = trues(1)
        # Only complete graphs with uniform edge weights have multiplicities [1, n - 1]
        elseif sort(last.(λ_counts)) == [1, n - 1]
            (
                eigvecs,
                cols_oneneg,
                bases_zerooneneg,
                bases_oneneg,
            ) = _completegraph_eigvecs_zerooneneg(λ_counts)
        else # All graphs except null, empty, and uniformly weighted complete graphs
            (
                eigvecs,
                cols_oneneg,
                bases_zerooneneg,
                bases_oneneg,
            ) = _generalgraph_eigvecs_zerooneneg(L, λ_counts)
        end
        
        return (eigvecs, cols_oneneg, bases_zerooneneg, bases_oneneg)
    end
    
    
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
    - `bases_zerooneneg::BitVector`: ADD LATER
    - `bases_oneneg::BitVector`: ADD LATER
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
        bases_zerooneneg = _column_span_basis.(eigvecs)
        bases_oneneg = _column_span_basis.(getindex.(eigvecs, :, cols_oneneg))
        
        # {-1,1}-bases are preferable to {-1,0,1}-bases whenever they exist
        idxs_replace = (length.(bases_zerooneneg) .== length.(bases_oneneg))
        bases_zerooneneg[idxs_replace] = copy(bases_oneneg[idxs_replace])
        
        return (eigvecs, cols_oneneg, bases_zerooneneg, bases_oneneg)
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
    - `bases_zerooneneg::BitVector`: ADD LATER
    - `bases_oneneg::BitVector`: ADD LATER
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
        
        # Test for linearly independent spanning {-1,0,1}- and {-1,1}-eigenbases
        bases_zerooneneg = _column_span_basis.(eigvecs)
        bases_oneneg = _column_span_basis.(getindex.(eigvecs, :, cols_oneneg))
        
        # {-1,1}-bases are preferable to {-1,0,1}-bases whenever they exist
        idxs_replace = (length.(bases_zerooneneg) .== length.(bases_oneneg))
        bases_zerooneneg[idxs_replace] = copy(bases_oneneg[idxs_replace])
        
        # bases_zerooneneg = trues(2) # Both eigenspaces have {-1,0,1}-bases for all n ≥ 2
        
        # if n % 2 == 0 # Both eigenspaces have {-1,1}-bases for all even n ≥ 2
        #     bases_oneneg = trues(2)
        # else # No odd-dimensional {-1,1}-vectors are orthogonal to the all-ones vector
        #     bases_oneneg = BitVector(undef, 2)
        #     bases_oneneg[idx_kernel] = true
        #     bases_oneneg[idx_nonkernel] = false
        # end
        
        return (eigvecs, cols_oneneg, bases_zerooneneg, bases_oneneg)
    end
    
    
    #- FUNCTION: `_get_cols_oneneg`
    """
        _get_cols_oneneg(eigvecs, idx_kernel, n)
    
    Find the column indices of all `{-1,1}`-eigenvectors of each matrix in `eigvecs`.
    
    Since every non-kernel `{-1,1}`-eigenvector must be orthogonal to the all-ones kernel
    vector and no odd-dimensional `{-1,1}`-vectors have an equal number of `-1`'s and `1`'s,
    no even-ordered Laplacian matriices are checked for non-kernel `{-1,1}`-eigenvectors.
    
    # Arguments
    - `eigvecs::Vector{Matrix{Int64}}`: an array of all `{-1,0,1}`-eigenvectors (as
        columns of matrices) in each eigenspace of some Laplacian matrix.
    - `idx_kernel::Int64`: the index of the `0`-eigenspace in `eigvecs`.
    - `n::Int64`: the order of the Laplacian whose eigenvectors `eigvecs` describes.
    
    # Returns
    - `cols_oneneg::Vector{Vector{Int64}}`: an array of all column indices of each
        matrix in `eigvecs` corresponding to `{-1,1}`-eigenvectors in that eigenspace.
    """
    function _get_cols_oneneg(eigvecs::Vector{Matrix{Int64}}, idx_kernel::Int64, n::Int64)
        # Only even-ordered Laplacians have non-kernel {-1,1}-eigenvectors
        if n % 2 == 0 # Check all eigenspaces for {-1,1}-eigenvectors
            cols_oneneg = findall.(v -> !any(v .== 0), eachcol.(eigvecs))
        else # Only check the kernel for {-1,1}-eigenvectors
            cols_oneneg = [Int64[] for _ in eigvecs]
            cols_oneneg[idx_kernel] = findall(
                v -> !any(v .== 0),
                eachcol(eigvecs[idx_kernel]),
            )
        end
        
        return cols_oneneg
    end
    
    
    #- FUNCTION: `_column_span_basis`
    """
        _column_span_basis(X)
    
    ADD LATER
    
    # Arguments
    - `X::Matrix{Int64}`: ADD LATER
    
    # Returns
    - `basis::Matrix{Int64}`: ADD LATER
    """
    function _column_span_basis(X::Matrix{Int64})
        reduced_form = rref_with_pivots(X)
        reduced_rank = rank(reduced_form[1], 1e-5)
        
        for comb in combinations(reduced_form[2], reduced_rank)
            if rank(X[:, comb], 1e-5) == reduced_rank
                basis = Int64.(X[:, comb])
                return basis
            end
        end
    end
end