module EigvecTests1
    #- EXPORTS AND IMPORTS
    export eigvecs_zerooneneg
    using LinearAlgebra: rank
    using Combinatorics: multiset_permutations
    
    
    #- FUNCTION: `eigvecs_zerooneneg`
    """
        eigvecs_zerooneneg(L::AbstractMatrix{Int64}, λ::Int64)
    
    Find all normalized {-1,0,1}- and {-1,1}-vectors in a Laplacian eigenspace.
    
    The first nonzero entry of each eigenvector is normalized to 1. If the graph order
    is even (or 1), identify the indices of all {-1,1}-vectors in the larger array.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of an undirected graph.
    - `λ::Int64`: the eigenvalue of some eigenspace of `L`.
    
    # Returns
    - `eigvecs::Vector{Vector{Int64}}`: an array of eigenvectors corresponding to `λ`.
    - `idxs_oneneg::Vector{Int64}`: the indices of all {-1,1}-eigenvectors in `eigvecs`.
    
    # Examples
    Test the 3-eigenspace of P_3 (not {-1,0,1}-diagonalizable):
    ```jldoctest
    julia> L = [ 1  -1   0;
                -1   2  -1;
                 0  -1   1];
    julia> eigvecs_zerooneneg(L, 3)
    (Vector{Int64}[], Int64[])
    ```
    
    Test the 3-eigenspace of K_3 (only {-1,0,1}-diagonalizable):
    ```jldoctest
    julia> L = [ 2  -1  -1;
                -1   2  -1;
                -1  -1   2];
    julia> eigvecs_zerooneneg(L, 3)
    ([[1, -1, 0], [1, 0, -1], [0, 1, -1]], Int64[])
    ```
    
    Test the 4-eigenspace of K_4 ({-1,1}-diagonalizable):
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
            else # `L` (which must have edge weights) has a higher nullity
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
        
        function fill_v!(v::Vector{Int64}, vec_slice2::Vector{Int64}, k::Int64)
            v[(k + 1):end] = vec_slice2
            return v
        end
        
        for k in 1:(n - 1) # Iterate over all possible indices of the first nonzero entry
            v = Vector{Int64}(undef, n)
            
            # Normalize the leading entry to 1
            v[1:(k - 1)] .= 0
            v[k] = 1
            
            r = n - k
            for i in 1:ceil(Int64, r/2) # Iterate over permutations of the remaining entries
                #= Every other eigenvector has an equal number of 1's and -1's (due to
                orthogonality with the all-ones vector of the kernel) =#
                entries_set = vcat(
                    -ones(Int64, i),
                    zeros(Int64, r - 2i + 1),
                    ones(Int64, i - 1)
                )
                
                # Place the current eigenvector generator in `generator_slices`
                generator_slices[idx] = (
                    fill_v!(v, vec_slice2, k)
                    for vec_slice2 in multiset_permutations(entries_set, r)
                )
                idx += 1
            end
        end
        
        return Iterators.flatten(generator_slices) # Chain all generators together
    end
end

module EigvecTests2
    #- EXPORTS AND IMPORTS
    export eigvecs_zerooneneg
    using LinearAlgebra: rank
    using Combinatorics: multiset_permutations
    
    
    #- FUNCTION: `eigvecs_zerooneneg`
    """
        eigvecs_zerooneneg(L::AbstractMatrix{Int64}, λ::Int64)
    
    Find all normalized {-1,0,1}- and {-1,1}-vectors in a Laplacian eigenspace.
    
    The first nonzero entry of each eigenvector is normalized to 1. If the graph order
    is even (or 1), identify the indices of all {-1,1}-vectors in the larger array.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of an undirected graph.
    - `λ::Int64`: the eigenvalue of some eigenspace of `L`.
    
    # Returns
    - `eigvecs::Vector{Vector{Int64}}`: an array of eigenvectors corresponding to `λ`.
    - `idxs_oneneg::Vector{Int64}`: the indices of all {-1,1}-eigenvectors in `eigvecs`.
    
    # Examples
    Test the 3-eigenspace of P_3 (not {-1,0,1}-diagonalizable):
    ```jldoctest
    julia> L = [ 1  -1   0;
                -1   2  -1;
                 0  -1   1];
    julia> eigvecs_zerooneneg(L, 3)
    (Vector{Int64}[], Int64[])
    ```
    
    Test the 3-eigenspace of K_3 (only {-1,0,1}-diagonalizable):
    ```jldoctest
    julia> L = [ 2  -1  -1;
                -1   2  -1;
                -1  -1   2];
    julia> eigvecs_zerooneneg(L, 3)
    ([[1, -1, 0], [1, 0, -1], [0, 1, -1]], Int64[])
    ```
    
    Test the 4-eigenspace of K_4 ({-1,1}-diagonalizable):
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
            else # `L` (which must have edge weights) has a higher nullity
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
            vec_slice1 = vcat(zeros(Int64, k - 1), 1) # Normalize the leading entry to 1
            r = n - k
            
            for i in 1:ceil(Int64, r/2) # Iterate over permutations of the remaining entries
                #= Every other eigenvector has an equal number of 1's and -1's (due to
                orthogonality with the all-ones vector of the kernel) =#
                entries_set = vcat(
                    -ones(Int64, i),
                    zeros(Int64, r - 2i + 1),
                    ones(Int64, i - 1)
                )
                
                # Place the current eigenvector generator in `generator_slices`
                generator_slices[idx] = (
                    vcat(vec_slice1, vec_slice2)
                    for vec_slice2 in multiset_permutations(entries_set, r)
                )
                idx += 1
            end
        end
        
        return Iterators.flatten(generator_slices) # Chain all generators together
    end
end

using BenchmarkTools

const L1 = [[14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 13  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1  0 13 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 11  0  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0 11  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0 11  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0  0 11 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1  8  0  0  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1  0  8  0  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1  0  0  8  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1  0  0  0  8  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  8  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  8  0]
    [-1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  0  8]];

const L2 = [[15 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 15 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 14  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1  0 14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 11  0  0  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0 11  0  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0 11  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0  0 11  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0  0  0 11 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  9  0  0  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  9  0  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  9  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  9  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  9  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  9  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  0  9]];

const K9 = [[ 8 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1  8 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1  8 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1  8 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  8 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1  8 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1  8 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1  8 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1  8]];

const K12 = [ 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    -1 11 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    -1 -1 11 -1 -1 -1 -1 -1 -1 -1 -1 -1;
    -1 -1 -1 11 -1 -1 -1 -1 -1 -1 -1 -1;
    -1 -1 -1 -1 11 -1 -1 -1 -1 -1 -1 -1;
    -1 -1 -1 -1 -1 11 -1 -1 -1 -1 -1 -1;
    -1 -1 -1 -1 -1 -1 11 -1 -1 -1 -1 -1;
    -1 -1 -1 -1 -1 -1 -1 11 -1 -1 -1 -1;
    -1 -1 -1 -1 -1 -1 -1 -1 11 -1 -1 -1;
    -1 -1 -1 -1 -1 -1 -1 -1 -1 11 -1 -1;
    -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 11 -1;
    -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 11];

const K15 = [[14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 14 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 14 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 14 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 14 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 14 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 14 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 14 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 14 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 14 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 14]];

const K18 = [[17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 17]];