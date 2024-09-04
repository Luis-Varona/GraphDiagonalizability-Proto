module ZeroOneNegBandwidths
    #- EXPORTS AND IMPORTS
    export is_DiagGraph
    
    using LinearAlgebra: rank
    using Combinatorics: combinations
    
    include("LaplacianSpectra.jl")
    include("KOrthogonalizability.jl")
    include("GraphObjects.jl")
    
    using .LaplacianSpectra: integer_eigvals_by_multi, eigvecs_zerooneneg
    using .KOrthogonalizability: is_k_orthogonalizable
    using .GraphObjects: DiagGraph, graph6_to_laplacian
    
    
    #- FUNCTION: `is_DiagGraph`
    """
        is_DiagGraph(L::AbstractMatrix{Int64}; min_zerooneneg = 1, min_oneneg = 1)
        is_DiagGraph(graph6_string::String; min_zerooneneg = 1, min_oneneg = 1)
    
    Determine whether an undirected graph is `{-1,0,1}`-diagonalizable.

    If `true`, also compute `{-1,0,1}`- and (if applicable) `{-1,1}`-eigendecompositions.
    (Data in stored in a `DiagGraph` object.) Optional lower bounds on the bandwidths are
    included as inputs in case the user already has theoretical bounds for a given graph.
    
    Weighted graphs are a supported format, but the edge weights must be integers.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of an undirected graph.
    - `graph6::String`: the graph6 string of the graph (if and only if `L` is not given).
    
    # Keywords
    - `min_zerooneneg::Int64 = 1`: the minimum `{-1,0,1}`-bandwidth for which to test. If
        only testing for diagonalizability and not the exact bandwidth, set this to `μ`.
    - `min_oneneg::Union{Float64, Int64} = 1`: the minimum `{-1,1}`-bandwidth for which to
        test. If only testing for diagonalizability and not the exact bandwidth, set this to
        `μ`. (If set to `Inf`, then `{-1,1}`-diagonalizability is not tested at all.)
    
    # Returns
    - `is_diag::Bool`: whether the graph is `{-1,0,1}`- (and thus `{-1,1}`-) diagonalizable.
    - `Γ::DiagGraph`: data on the `{-1,0,1}`- and `{-1,1}`-spectra of the graph.
    
    # Examples
    Test the (unweighted) complete graph `K_10` in graph6 format (with theoretical lower
        bounds on the S-bandwidths of complete graphs):
    ```jldoctest
    julia> s = "I~~~~~~~w";
    julia> (is_diag, Γ) = is_DiagGraph(s, min_zerooneneg = 2, min_oneneg = 9);
    julia> (is_diag, Γ.band_zerooneneg, Γ.band_oneneg)
    (true, {-1,0,1}-diagonalizable graph on 10 vertices, 2, 9)
    julia> Γ.eigvals
    10-element Vector{Int64}:
      0
     10
     10
     10
     10
     10
     10
     10
     10
     10
    julia> Γ.eigvecs_zerooneneg
    10×10 Matrix{Int64}:
     1   0   1   1   1   1   1   0   0   0
     1   0   1   1   0  -1  -1   0   0   0
     1   0   1   0  -1   0   1   0   0   0
     1   0   1   0   0   0  -1   1   0   0
     1   0  -1  -1   0   0   0   0   1   1
     1   0  -1  -1   0   0   0   0  -1  -1
     1   0  -1   0   0   0   0   0  -1   0
     1   0  -1   0   0   0   0   1   1   0
     1   1   0   0   0   0   0  -1   0   0
     1  -1   0   0   0   0   0  -1   0   0
    julia> Γ.eigvecs_oneneg
    10×10 Matrix{Int64}:
     1   1   1   1   1   1   1   1   1   1
     1  -1  -1  -1  -1  -1  -1  -1  -1   1
     1  -1  -1  -1  -1  -1  -1  -1   1  -1
     1  -1  -1  -1  -1  -1  -1   1  -1  -1
     1  -1  -1  -1  -1  -1   1  -1  -1  -1
     1  -1   1   1   1   1  -1  -1  -1  -1
     1   1  -1   1   1   1  -1  -1  -1  -1
     1   1   1  -1   1   1   1   1   1   1
     1   1   1   1  -1   1   1   1   1   1
     1   1   1   1   1  -1   1   1   1   1
    ```
    
    Test the (unweighted) complete multipartite graph `K_{1,2,3,5}` in Laplacian format:
    ```jldoctest
    julia> L = [ 6   0   0   0   0  -1  -1  -1  -1  -1  -1;
                 0   6   0   0   0  -1  -1  -1  -1  -1  -1;
                 0   0   6   0   0  -1  -1  -1  -1  -1  -1;
                 0   0   0   6   0  -1  -1  -1  -1  -1  -1;
                 0   0   0   0   6  -1  -1  -1  -1  -1  -1;
                -1  -1  -1  -1  -1   8   0   0  -1  -1  -1;
                -1  -1  -1  -1  -1   0   8   0  -1  -1  -1;
                -1  -1  -1  -1  -1   0   0   8  -1  -1  -1;
                -1  -1  -1  -1  -1  -1  -1  -1   9   0  -1;
                -1  -1  -1  -1  -1  -1  -1  -1   0   9  -1;
                -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  10];
    julia> (is_diag, Γ) = is_DiagGraph(L);
    julia> (is_diag, Γ, Γ.band_zerooneneg, Γ.band_oneneg)
    (false, {-1,0,1}-non-diagonalizable graph on 11 vertices, Inf, Inf)
    julia> Γ.eigvals
    11-element Vector{Int64}:
      0
      9
      8
      8
     11
     11
     11
      6
      6
      6
      6
    julia> Γ.eigvecs_zerooneneg
    11×0 Matrix{Int64}
    julia> Γ.eigvecs_oneneg
    11×0 Matrix{Int64}
    ```
    
    Test a weighted version of the strong product `K_{1,1,2} ⊠ K_{2,2}` in Laplacian format
        (here the spectra are not shown for the sake of brevity):
    ```jldoctest
    julia> using LinearAlgebra: kron
    julia> K_112 = [ 2   0  -1  -1;
                     0   2  -1  -1;
                    -1  -1   3  -1;
                    -1  -1  -1   3];
    julia> K_22 = [ 2  -1  -1   0;
                   -1   2   0  -1;
                   -1   0   2  -1;
                    0  -1  -1   2];
    julia> L = kron(K_112, K_22);
    julia> (is_diag, Γ) = is_DiagGraph(L);
    julia> (is_diag, Γ, Γ.band_zerooneneg, Γ.band_oneneg)
    (true, {-1,0,1}-diagonalizable graph on 16 vertices, 1, Inf)
    ```
    """
    function is_DiagGraph(
        graph6::String;
        kwargs...
    )::Tuple{Bool, DiagGraph}
        
        L = graph6_to_laplacian(graph6)
        return is_DiagGraph(L; kwargs...)
    end
    
    
    function is_DiagGraph(
        L::AbstractMatrix{Int64};
        min_zerooneneg::Int64 = 1,
        min_oneneg::Union{Float64, Int64} = 1
    )::Tuple{Bool, DiagGraph}
        
        n = size(L, 2)
        (is_integral, Λ, M) = integer_eigvals_by_multi(L)
        
        # Assuming integer weights, {-1,0,1}-diagonalizability implies Laplacian integrality
        if !is_integral
            is_diag = false
            Γ = DiagGraph(
                L, Inf, Inf, Λ, Matrix{Int64}(undef, n, 0), Matrix{Int64}(undef, n, 0)
            )
        else
            V_zerooneneg = Matrix{Int64}(undef, n, n)
            V_oneneg = copy(V_zerooneneg) # Use a new pointer
            j = 1
            
            for (λ, μ) in M
                begin # Only test for bandwidths no less than those of previous eigenspaces
                    (min_zerooneneg, min_oneneg, basis_zerooneneg, basis_oneneg) =
                    eigspace_bandwidths(
                        L, λ, μ, min_zerooneneg = min_zerooneneg, min_oneneg = min_oneneg
                    )
                end
                
                # If any eigenspace has infinite bandwidth, the graph is not diagonalizable
                if min_zerooneneg == Inf
                    is_diag = false
                    Γ = DiagGraph(
                        L, Inf, Inf, Λ,
                        Matrix{Int64}(undef, n, 0), Matrix{Int64}(undef, n, 0)
                    )
                    break
                end
                
                # Combine the computed eigenbases to (hopefully) form diagonalizing matrices
                V_zerooneneg[:, j:(j + μ - 1)] = basis_zerooneneg
                if min_oneneg != Inf
                    V_oneneg[:, j:(j + μ - 1)] = basis_oneneg
                end
                
                j += μ
            end
            
            if min_zerooneneg != Inf
                if min_oneneg == Inf
                    V_oneneg = Matrix{Int64}(undef, n, 0)
                end
                
                # All checks have been passed, so the graph is {-1,0,1}-diagonalizable
                is_diag = true
                Γ = DiagGraph(L, min_zerooneneg, min_oneneg, Λ, V_zerooneneg, V_oneneg)
            end
        end
        
        return (is_diag, Γ)
    end
    
    
    #- FUNCTION: `eigspace_bandwidths`
    """
        eigspace_bandwidths(L, λ, μ; min_zerooneneg = 1, min_oneneg = 1)
    
    Compute the `{-1,0,1}`- and `{-1,1}`-bandwidths of a Laplacian eigenspace.
    
    Optional lower bounds on the bandwidths are included to avoid testing any
    eigenspace for a bandwidth less than the minimum of a previous eigenspace,
    and in case the user already has theoretical bounds for a given graph.
    
    Weighted Laplacians are a supported format, but the edge weights must be integers.
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: the Laplacian matrix of an undirected graph.
    - `λ::Int64`: the eigenvalue of some eigenspace of `L`.
    - `μ::Int64`: the dimension of the eigenspace.
    
    # Keywords
    - `min_zerooneneg::Int64 = 1`: the minimum `{-1,0,1}`-bandwidth for which to test.
    - `min_oneneg::Union{Float64, Int64} = 1`: the minimum `{-1,1}`-bandwidth for which to
        test. (If set to `Inf`, then the `{-1,1}`-bandwidth is not tested.)
    
    # Returns
    - `band_zerooneneg::Union{Float64, Int64}`: the `{-1,0,1}`-bandwidth of the eigenspace.
    - `band_oneneg::Union{Float64, Int64}`: the `{-1,1}`-bandwidth of the eigenspace.
    - `basis_zerooneneg::Matrix{Int64}`: a `{-1,0,1}`-basis with minimized bandwidth.
    - `basis_oneneg::Matrix{Int64}`: a `{-1,1}`-basis with minimized bandwidth.
    """
    function eigspace_bandwidths(
        L::AbstractMatrix{Int64},
        λ::Int64,
        μ::Int64;
        min_zerooneneg::Int64 = 1,
        min_oneneg::Union{Float64, Int64} = 1
    )::Tuple{Union{Float64, Int64}, Union{Float64, Int64}, Matrix{Int64}, Matrix{Int64}}
        
        n = size(L, 1)
        (eigvecs, idxs_oneneg) = eigvecs_zerooneneg(L, λ)
        num_eigvecs = length(eigvecs)
        
        if num_eigvecs < μ # There are too few vectors to form an {-1,0,1}-basis
            band_zerooneneg = Inf
            band_oneneg = Inf # Always greater than or equal to the {-1,0,1}-bandwidth
            basis_zerooneneg = Matrix{Int64}(undef, n, 0)
            basis_oneneg = copy(basis_zerooneneg) # Use a new pointer
        else
            if length(idxs_oneneg) < μ # There are too few vectors to form an {-1,1}-basis
                # Set this instead of `band_oneneg` to control later if-else statements
                min_oneneg = Inf
            end
            
            if μ == 1 # All m × 1 matrices are pairwise orthogonal
                band_zerooneneg = min_zerooneneg
                
                if min_oneneg == Inf
                    band_oneneg = Inf
                    basis_zerooneneg = reshape(eigvecs[1], n, 1)
                    basis_oneneg = Matrix{Int64}(undef, n, 0)
                else
                    band_oneneg = max(min_oneneg, band_zerooneneg)
                    basis_zerooneneg = reshape(eigvecs[idxs_oneneg[1]], n, 1)
                    basis_oneneg = copy(basis_zerooneneg) # Use a new pointer
                end
            else
                """
                    DFS(idxs, k, idx_set)
                
                Recursively search for a complete `{-1,0,1}`- or `{-1,1}`-eigenbasis.
                
                Start with all 2-combinations of an array of eigenvectors and recursively
                add subsequent vectors to find a k-orthogonal basis set. (Array indices are
                used instead of the vectors themselves to conserve memory.)
                
                # Arguments
                - `idxs::Vector{Int64}`: the current set of indices.
                - `k::Int64`: the `{-1,0,1}`- or `{-1,1}`-bandwidth for which to test.
                - `idx_set::Vector{Int64}`: the index set (all eigenvectors if searching for
                    a `{-1,0,1}`-basis and only a subset if searching for a `{-1,1}`-basis).
                
                # Returns
                - `is_k_ortho::Bool`: whether a `k`-orthogonal basis set was found.
                - `idxs_final::Vector{Int64}`: the indices of the `k`-orthogonal basis set.
                """
                function DFS(
                    idxs::Vector{Int64},
                    k::Int64,
                    idx_set::Vector{Int64}
                )::Tuple{Bool, Vector{Int64}}
                    
                    depth = length(idxs)
                    eigbasis = hcat(eigvecs[idxs]...) # The current eigenbasis submatrix
                    
                    if rank(eigbasis, 1e-5) != depth # Verify linear independence
                        is_k_ortho = false
                        idxs_final = Int64[]
                    else
                        (is_k_ortho, order) = is_k_orthogonalizable(eigbasis, k)
                        
                        if !is_k_ortho
                            idxs_final = Int64[]
                        elseif depth == μ # A complete set of basis vectors has been found
                            idxs_final = idxs[order] # Impose a k-orthogonal column ordering
                        else
                            is_k_ortho = false
                            idxs_final = Int64[]
                            last = idxs[end] # The greatest index in the (sorted) list
                            
                            # Only add subsequent indices to avoid duplicate combinations
                            for idx_next in filter(x -> x > last, idx_set)
                                idxs_next = vcat(idxs, idx_next)
                                
                                # Recursively call `DFS` and search for a complete basis set
                                (is_k_ortho, idxs_final) = DFS(idxs_next, k, idx_set)
                                if is_k_ortho
                                    break
                                end
                            end
                        end
                    end
                    
                    return (is_k_ortho, idxs_final)
                end
                
                """
                    search_roots(roots, min_k, idx_set)
                
                Iteratively call `DFS` to find the exact `{-1,0,1}`- or `{-1,1}`-bandwidth.
                
                # Arguments
                - `idx_set::Vector{Int64}`: the index set of all eigenvectors to search.
                - `min_k::Int64`: the minimum eigenspace bandwidth to test.
                
                # Returns
                - `k::Union{Float64, Int64}`: the bandwidth of the eigenspace.
                - `basis::Matrix{Int64}`: an eigenbasis with minimized bandwidth.
                """
                function search_roots(
                    idx_set::Vector{Int64},
                    min_k::Int64
                )::Tuple{Union{Float64, Int64}, Matrix{Int64}}
                    
                    has_basis = false
                    idxs_temp = Int64[]
                    roots = combinations(idx_set, 2) # All 2-combinations of eigenvectors
                    
                    # Verify that a μ-orthogonal basis exists at all before continuing
                    for root in roots
                        # Every m × μ eigenbasis is automatically μ-orthogonal
                        (has_basis, idxs_temp) = DFS(root, μ, idx_set)
                        if has_basis
                            break
                        end
                    end
                    
                    if !has_basis # Lower bandwidths are even stronger conditions
                        k = Inf
                        basis = Matrix{Int64}(undef, n, 0)
                    else
                        has_basis = false
                        idxs_basis = Int64[]
                        k = min_k - 1 # Only begin testing at the lower bandwidth bound
                        
                        # Iteratively search for a (k + 1)-eigenbasis until one is found
                        while !has_basis && k < max(μ, min_k) - 1
                            k += 1
                            for root in roots
                                # Call `DFS` on some some 2-combination of eigenvectors
                                (has_basis, idxs_basis) = DFS(root, k, idx_set)
                                if has_basis
                                    break
                                end
                            end
                        end
                        
                        if has_basis
                            basis = hcat(eigvecs[idxs_basis]...) 
                        else # Use the initially computed μ-orthogonal eigenbasis
                            k = max(μ, min_k)
                            basis = hcat(eigvecs[idxs_temp]...)
                        end
                    end
                    
                    return (k, basis)
                end
                
                # Compute the {-1,0,1}-bandwidth by searching all {-1,0,1}-eigenvectors
                (band_zerooneneg, basis_zerooneneg) = search_roots(
                    collect(1:num_eigvecs), min_zerooneneg
                )
                
                if min_oneneg == Inf
                    band_oneneg = Inf
                    basis_oneneg = Matrix{Int64}(undef, n, 0)
                else
                    if band_zerooneneg == Inf
                        band_oneneg = Inf
                        basis_oneneg = Matrix{Int64}(undef, n, 0)
                    elseif !any(basis_zerooneneg .== 0) # Also a valid {-1,1}-eigenbasis
                        band_oneneg = max(min_oneneg, band_zerooneneg)
                        basis_oneneg = copy(basis_zerooneneg) # Use a new pointer
                    else
                        # Compute the {-1,1}-bandwidth by searching only {-1,1}-eigenvectors
                        (band_oneneg, basis_oneneg) = search_roots(idxs_oneneg, min_oneneg)
                        if band_oneneg == band_zerooneneg # The {-1,1}-basis is preferable
                            basis_zerooneneg = copy(basis_oneneg) # Use a new pointer
                        end
                    end
                end
            end
        end
        
        return (band_zerooneneg, band_oneneg, basis_zerooneneg, basis_oneneg)
    end
end