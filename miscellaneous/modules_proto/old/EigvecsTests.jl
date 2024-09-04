#- #- #-
module A
    #- EXPORTS AND IMPORTS
    export eigvecs_zerooneneg
    using LinearAlgebra: rank
    using Combinatorics: multiset_permutations
    
    #- FUNCTION: `eigvecs_zerooneneg`
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
    function laplacian_eigvecs(n::Int64)::Iterators.Flatten{Vector{Base.Generator}}
        generator_slices = Vector{Base.Generator}(undef, floor(Int64, n^2/4))
        idx = 1
        
        for k in 1:(n - 1) # Iterate over all possible indices of the first nonzero entry
            leading = vcat(zeros(Int64, k - 1), 1) # Normalize the leading entry to 1
            r = n - k
            
            for i in 1:ceil(Int64, r/2) # Iterate over permutations of the remaining entries
                #= Every other eigenvector has an equal number of 1's and -1's (due to
                orthogonality with the all-ones vector of the kernel) =#
                entries_set = vcat(
                    -ones(Int64, i), zeros(Int64, r - 2i + 1), ones(Int64, i - 1)
                )
                
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



#- #- #-
module B
    #- EXPORTS AND IMPORTS
    export eigvecs_zerooneneg
    using LinearAlgebra: rank
    using Combinatorics: multiset_permutations
    
    #- FUNCTION: `eigvecs_zerooneneg`
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
    function laplacian_eigvecs(n::Int64)::Iterators.Flatten{Vector{Base.Generator}}
        generator_slices = Vector{Base.Generator}(undef, floor(Int64, n^2/4))
        idx = 1
        
        for k in 1:(n - 1) # Iterate over all possible indices of the first nonzero entry
            leading = vcat(zeros(Int64, k - 1), 1) # Normalize the leading entry to 1
            r = n - k
            
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



#- #- #-
module C
    #- EXPORTS AND IMPORTS
    export eigvecs_zerooneneg
    
    using StaticArrays
    using LinearAlgebra: rank
    using Combinatorics: multiset_permutations
    
    #- FUNCTION: `eigvecs_zerooneneg`
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
    function laplacian_eigvecs(n::Int64)::Iterators.Flatten{Vector{Base.Generator}}
        generator_slices = Vector{Base.Generator}(undef, floor(Int64, n^2/4))
        idx = 1
        
        for k in 1:(n - 1) # Iterate over all possible indices of the first nonzero entry
            leading = vcat(zeros(Int64, k - 1), 1) # Normalize the leading entry to 1
            r = n - k
            
            for i in 1:ceil(Int64, r/2) # Iterate over permutations of the remaining entries
                #= Every other eigenvector has an equal number of 1's and -1's (due to
                orthogonality with the all-ones vector of the kernel) =#
                entries_set = SV(vcat(
                    -ones(Int64, i), zeros(Int64, r - 2i + 1), ones(Int64, i - 1)
                ))
                
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



#- #- #-
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

const K10_minus = [[ 8  0 -1 -1 -1 -1 -1 -1 -1 -1]
    [ 0  8 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1  9 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1  9 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  9 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1  9 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1  9 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1  9 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1  9 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  9]];

const K11_minus2 = [[ 9 -1 -1  0 -1 -1 -1 -1 -1 -1 -1]
    [-1  9  0 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1  0  9 -1 -1 -1 -1 -1 -1 -1 -1]
    [ 0 -1 -1  9 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 10 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 10 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 10 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 10 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 10 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 10 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 10]];