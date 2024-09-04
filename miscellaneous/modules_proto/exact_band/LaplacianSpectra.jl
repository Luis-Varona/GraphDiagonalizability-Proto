module LaplacianSpectra
    #- EXPORTS AND IMPORTS
    export integer_eigvals_by_multi, eigvecs_zerooneneg
    
    using LinearAlgebra: Symmetric, eigvals, rank
    using Combinatorics: multiset_permutations
    using DataStructures: counter
    
    
    #- proto eigvecs with loading
    function eigvecs_zerooneneg(
        L::AbstractMatrix{Int64},
        λ::Int64,
        E::Matrix{Int64}
    )::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
        
        if λ == 0
            n = size(L, 2)
            if rank(L, 1e-5) == n - 1 # The kernel is fully spanned by the all-ones vector
                eigvecs = [ones(Int64, n)]
                idxs_oneneg = [1]
            else # `L` (which must be weighted and/or disconnected) has a larger kernel
                (eigvecs, idxs_oneneg) = get_norm_eigvecs(L, λ, n, E)
                
                # Manually add the all-ones vector and update the {-1,1}-indices
                eigvecs = [ones(Int64, n), eigvecs...]
                idxs_oneneg = [1, (idxs_oneneg .+ 1)...]
            end
        else # The general case
            (eigvecs, idxs_oneneg) = get_norm_eigvecs(L, λ, size(L, 2), E)
        end
        
        return (eigvecs, idxs_oneneg)
    end
    
    function get_norm_eigvecs(
        L::AbstractMatrix{Int64},
        λ::Int64,
        n::Int64,
        E::Matrix{Int64}
    )::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
        
        eigvecs = Vector{Int64}[]
        idxs_oneneg = Int64[]
        
        if n % 2 == 0 # Only K_1 and even-ordered graphs are {-1,1}-diagonalizable
            for v in E # Generate possible eigenvectors
                if L*v == λ*v
                    push!(eigvecs, v)
                    if !any(v .== 0) # Check which are {-1,1}-eigenvectors
                        push!(idxs_oneneg, length(eigvecs))
                    end
                end
            end
        else # Do not check odd-ordered graphs for {-1,1}-eigenvectors
            for v in E
                if L*v == λ*v
                    push!(eigvecs, v)
                end
            end
        end
        
        return (eigvecs, idxs_oneneg)
    end
    
    
    #- FUNCTION: `integer_eigvals_by_multi`
    """
        integer_eigvals_by_multi(L::AbstractMatrix{Int64})
    
    Determine whether an undirected graph is Laplacian integral and find its spectra if so.
    
    Weighted Laplacians are a supported format, but the edge weights must be integers.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of an undirected graph.
    
    # Returns
    - `is_integral::Bool`: whether the eigenvalues are all integers.
    - `Λ::Vector{Int64}`: if `is_integral`, an array of all eigenvalues sorted by ascending
        multiplicity. Otherwise, an empty array.
    - `M::Vector{Pair{Int64, Int64}}`: if `is_integral`, an array of unique eigenvalues and
        their multiplicities. Otherwise, an empty array.
    
    # Examples
    Test the (unweighted) cycle graph `C_6`:
    ```jldoctest
    julia> L = [ 2  -1  -1   0   0   0;
                -1   2   0   0  -1   0;
                -1   0   2   0   0  -1;
                 0   0   0   2  -1  -1;
                 0  -1   0  -1   2   0;
                 0   0  -1  -1   0   2];
    julia> integer_eigvals_by_multi(L)
    (true, [0, 4, 3, 3, 1, 1], [0 => 1, 4 => 1, 3 => 2, 1 => 2])
    ```
    
    Test a weighted version of the star graph `S_5`:
    ```jldoctest
    julia> L = [ 1   0   0   0  -1;
                 0   1   0   0  -1;
                 0   0   1   0  -1;
                 0   0   0   7  -7;
                -1  -1  -1  -7  10];
    julia> integer_eigvals_by_multi(L)
    (false, Int64[], Pair{Int64, Int64}[])
    ```
    
    Test a weighted version of the strong product `K_2 ⊠ K_{1,1,2}`:
    ```jldoctest
    julia> L = [ 6  -2   0   0  -1  -1  -1  -1;
                -2   6   0   0  -1  -1  -1  -1;
                 0   0   6  -2  -1  -1  -1  -1;
                 0   0  -2   6  -1  -1  -1  -1;
                -1  -1  -1  -1   9  -3  -1  -1;
                -1  -1  -1  -1  -3   9  -1  -1;
                -1  -1  -1  -1  -1  -1   9  -3;
                -1  -1  -1  -1  -1  -1  -3   9];
    julia> integer_eigvals_by_multi(L)
    (true, [0, 4, 12, 12, 8, 8, 8, 8], [0 => 1, 4 => 1, 12 => 2, 8 => 4])
    ```
    """
    function integer_eigvals_by_multi(
        L::AbstractMatrix{Int64}
    )::Tuple{Bool, Vector{Int64}, Vector{Pair{Int64, Int64}}}
        
        # Computing eigenvalues is more efficient for `Symmetric`-type matrices
        Λ_float = eigvals(Symmetric(L, :L))
        Λ = round.(real.(Λ_float))
        is_integral = isapprox(Λ, Λ_float)
        
        if !is_integral
            Λ = Int64[]
            λ_counts = Pair{Int64, Int64}[]
        else
            # Count all multiplicitiess and sort unique eigenvalues by multiplicity
            λ_counts = sort(collect(counter(Int64.(Λ))), by = last)
            
            # Sort all eigenvalues, including repetitions, by multiplicity
            Λ = vcat(fill.(first.(λ_counts), last.(λ_counts))...)
        end
        
        return (is_integral, Λ, λ_counts)
    end
    
    
    #- FUNCTION: `eigvecs_zerooneneg`
    """
        eigvecs_zerooneneg(L::AbstractMatrix{Int64}, λ::Int64)
    
    Find all normalized {-1,0,1}- and {-1,1}-vectors in a Laplacian eigenspace.
    
    Weighted Laplacians are a supported format, but the edge weights must be integers.
    
    The first nonzero entry of each eigenvector is normalized to 1. If the graph order
    is even (or 1), identify the indices of all {-1,1}-vectors in the larger array.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of an undirected graph.
    - `λ::Int64`: the eigenvalue of some eigenspace of `L`.
    
    # Returns
    - `eigvecs::Vector{Vector{Int64}}`: an array of eigenvectors corresponding to `λ`.
    - `idxs_oneneg::Vector{Int64}`: the indices of all {-1,1}-eigenvectors in `eigvecs`.
    
    # Examples
    Test the `3`-eigenspace of the (unweighted) path graph `P_3`:
    ```jldoctest
    julia> L = [ 1  -1   0;
                -1   2  -1;
                 0  -1   1];
    julia> eigvecs_zerooneneg(L, 3)
    (Vector{Int64}[], Int64[])
    ```
    
    Test the `4`-eigenspace of a weighted version of the strong product `P_2 ⊠ P_3`:
    ```jldoctest
    julia> L = [ 3  -1   0   0  -1  -1;
                -1   3   0   0  -1  -1;
                 0   0   3  -1  -1  -1;
                 0   0  -1   3  -1  -1;
                -1  -1  -1  -1   6  -2;
                -1  -1  -1  -1  -2   6];
    julia> eigvecs_zerooneneg(L, 4)
    ([[1, -1, 0, 0, 0, 0], [1, -1, -1, 1, 0, 0], [1, -1, 1, -1, 0, 0], [0, 0, 1, -1, 0, 0]], Int64[])
    ```
    
    Test the `4`-eigenspace of the (unweighted) complete graph `K_4`:
    ```jldoctest
    julia> L = [ 3  -1  -1  -1;
                -1   3  -1  -1;
                -1  -1   3  -1;
                -1  -1  -1   3];
    julia> eigvecs_zerooneneg(L, 4)
    ([[-1, 1, 1, 1], [1, -1, 1, 1], [1, 1, -1, 1], [1, 1, 1, -1]], [1, 2, 3, 4])
    ```
    """
    function eigvecs_zerooneneg(
        L::AbstractMatrix{Int64},
        λ::Int64
    )::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
        
        n = size(L, 2)
        
        if n == 1 # The sole eigenspace of K_1 is spanned by the all-ones vector
            eigvecs = [ones(Int64, n)]
            idxs_oneneg = [1]
        elseif λ == 0
            if rank(L, 1e-5) == n - 1 # The kernel is fully spanned by the all-ones vector
                eigvecs = [ones(Int64, n)]
                idxs_oneneg = [1]
            else # `L` (which must be weighted and/or disconnected) has a larger kernel
                (eigvecs, idxs_oneneg) = get_norm_eigvecs(L, λ, n)
                
                # Manually add the all-ones vector and update the {-1,1}-indices
                eigvecs = [ones(Int64, n), eigvecs...]
                idxs_oneneg = [1, (idxs_oneneg .+ 1)...]
            end
        else # The general case
            (eigvecs, idxs_oneneg) = get_norm_eigvecs(L, λ, n)
        end
        
        return (eigvecs, idxs_oneneg)
    end
    
    
    #- FUNCTION: `get_norm_eigvecs`
    """
    get_norm_eigvecs()
    
    Find all normalized {-1,0,1}-eigenvectors except the all-ones vector, if applicable.
    
    If the graph order is even (or 1), identify the indices of all {-1,1}-vectors as well.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of an undirected graph.
    - `λ::Int64`: the eigenvalue of some eigenspace of `L`.
    - `n::Int64`: the order of the undirected graph represented by `L`. (Since
        `get_norm_eigvecs` is only called by `eigvecs_zerooneneg`, which already computes
        `n`, passing this argument avoids redundant computation.)
    
    # Returns
    - `eigvecs::Vector{Vector{Int64}}`: an array of eigenvectors corresponding to `λ`.
    - `idxs_oneneg::Vector{Int64}`: the indices of all {-1,1}-eigenvectors in `eigvecs`.
    """
    function get_norm_eigvecs(
        L::AbstractMatrix{Int64},
        λ::Int64,
        n::Int64,
    )::Tuple{Vector{Vector{Int64}}, Vector{Int64}}
        
        eigvecs = Vector{Int64}[]
        idxs_oneneg = Int64[]
        
        if n % 2 == 0 # Only K_1 and even-ordered graphs are {-1,1}-diagonalizable
            for v in laplacian_eigvecs(n) # Generate possible eigenvectors
                if L*v == λ*v
                    push!(eigvecs, v)
                    if !any(v .== 0) # Check which are {-1,1}-eigenvectors
                        push!(idxs_oneneg, length(eigvecs))
                    end
                end
            end
        else # Do not check odd-ordered graphs for {-1,1}-eigenvectors
            for v in laplacian_eigvecs(n)
                if L*v == λ*v
                    push!(eigvecs, v)
                end
            end
        end
        
        return (eigvecs, idxs_oneneg)
    end
    
    
    #- FUNCTION: `laplacian_eigvecs`
    """
        laplacian_eigvecs(n::Int64)
    
    Return a lazy generator of all n-dimensional Laplacian {-1,0,1}-eigenvectors.
    
    Each eigenvector is normalized such that its first nonzero entry is 1 (thus enforcing
    pairwise independence between all generated vectors).
    
    # Arguments
    - `n::Int64`: the order of the undirected graph represented by some Laplacian matrix.
    
    # Returns
    - `Iterators.Flatten{Vector{Base.Generator}}`: a lazily evaluated iterator of all
        possible {-1,0,1}-eigenvectors of some undirected Laplacian of order `n`.
    """
    function laplacian_eigvecs(n::Int64)::Iterators.Flatten{Vector{Base.Generator}}
        generator_slices = Vector{Base.Generator}(undef, floor(Int64, n^2/4))
        idx = 1
        
        for k in 1:(n - 1) # Iterate over all possible indices of the first nonzero entry
            # Normalize the leading (first nonzero) entry to 1
            leading = Vector{Int64}(undef, k)
            leading[1:(k - 1)] .= 0
            leading[k] = 1
            
            r = n - k # The number of remaining entries in the eigenvector
            for i in 1:ceil(Int64, r/2) # Iterate over permutations of the remaining entries
                #= Every other eigenvector has an equal number of 1's and -1's (due to
                orthogonality with the all-ones vector of the kernel) =#
                entries_set = Vector{Int64}(undef, r)
                entries_set[1:i] .= -1
                entries_set[(i + 1):(r - i + 1)] .= 0
                entries_set[(r - i + 2):r] .= 1
                
                # Place the current eigenvector generator in `generator_slices`
                generator_slices[idx] = (
                    vcat(leading, entries)
                    for entries in multiset_permutations(entries_set, r)
                )
                idx += 1
            end
        end
        
        return Iterators.flatten(generator_slices) # Chain all generators together
    end
end