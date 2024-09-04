module EigenvectorGenerators
    # EXPORTS AND IMPORTS
    export potential_kernel_eigvecs
    export potential_nonkernel_eigvecs
    using Combinatorics: multiset_permutations
    
    
    #- FUNCTION: `potential_kernel_eigvecs`
    """
        potential_kernel_eigvecs(n)
    
    Return a lazy generator of potential kernel `{-1,0,1}`-eigenvectors of a Laplacian.
    
    Each vector is normalized so that its first nonzero entry is `1`, enforcing pairwise
    linear independence between all generated vectors.
    
    # Arguments
    - `n::Int64`: the order of the Laplacian matrix of some undirected graph for which to
        find potential kernel vectors.
    
    # Returns
    - `eigvec_generator::Iterators.Flatten{Vector{Base.Generator}}`: a lazily evaluated
        iterator over all `n`-dimensional `{-1,0,1}`-vectors, unique up to span.
    
    # Examples
    Generate all potential kernel eigenvectors for an order `3` Laplacian matrix:
    ```jldoctest
    julia> collect(Vector{Int64}, potential_kernel_eigvecs(3))
    13-element Vector{Vector{Int64}}:
     [1, 1, 1]
     [1, -1, 1]
     [1, 0, 1]
     [1, 1, -1]
     [1, -1, -1]
     [1, 0, -1]
     [1, 1, 0]
     [1, -1, 0]
     [1, 0, 0]
     [0, 1, 1]
     [0, 1, -1]
     [0, 1, 0]
     [0, 0, 1]
    ```
    """
    function potential_kernel_eigvecs(n::Int64)
        generator_slices = Vector{Base.Generator}(undef, n)
        
        for k in 1:n # Iterate over all possible indices of the first nonzero entry
            leading = Vector{Int64}(undef, k)
            leading[1:(k - 1)] .= 0
            leading[k] = 1 # Normalize this leading entry to 1
            
            # The remaining entries are (n - k)-permutations with repetition of {-1,0,1}
            entries_perms = Iterators.product(Iterators.repeated((1, -1, 0), n - k)...)
            # Add a generator of all such possible permutations to `generator_slices`
            generator_slices[k] = (vcat(leading, body...) for body in entries_perms)
        end
        
        # Chain all generators together to form a single lazy iterator
        eigvec_generator = Iterators.flatten(generator_slices)
        return eigvec_generator
    end
    
    
    #- FUNCTION: `potential_nonkernel_eigvecs`
    """
        potential_nonkernel_eigvecs(n)
    
    Return a lazy generator of potential non-kernel `{-1,0,1}`-eigenvectors of a Laplacian.
    
    Each vector is normalized so that its first nonzero entry is `1`, enforcing pairwise
    linear independence between all generated vectors. Since all Laplacian matrices have
    pairwise orthogonal eigenspaces and the all-ones vector is always in the kernel, every
    non-kernel `{-1,0,1}`-eigenvector must also have an equal number of `-1`'s and `1`'s.
    
    # Arguments
    - `n::Int64`: the order of the Laplacian matrix of some undirected graph for which to
        find potential non-kernel eigenvectors.
    
    # Returns
    - `eigvec_generator::Iterators.Flatten{Vector{Base.Generator}}`: a lazily evaluated
        iterator over all `n`-dimensional `{-1,0,1}`-vectors orthogonal to the all-ones
        kernel vector, unique up to span.
    
    # Examples
    Generate all potential non-kernel eigenvectors of an order `4` Laplacian matrix:
    ```jldoctest
    julia> collect(Vector{Int64}, potential_nonkernel_eigvecs(4))
    9-element Vector{Vector{Int64}}:
     [1, -1, 0, 0]
     [1, 0, -1, 0]
     [1, 0, 0, -1]
     [1, -1, -1, 1]
     [1, -1, 1, -1]
     [1, 1, -1, -1]
     [0, 1, -1, 0]
     [0, 1, 0, -1]
     [0, 0, 1, -1]
    ```
    """
    function potential_nonkernel_eigvecs(n::Int64)
        generator_slices = Vector{Base.Generator}(undef, floor(Int64, n^2/4))
        idx = 1
        
        for k in 1:(n - 1) # Iterate over all possible indices of the first nonzero entry
            leading = Vector{Int64}(undef, k)
            leading[1:(k - 1)] .= 0
            leading[k] = 1 # Normalize this leading entry to 1
            r = n - k # The number of remaining entries in the eigenvector
            
            # Iterate over permutations of remaining entries, varying the number of -1's
            for i in 1:ceil(Int64, r/2)
                entries_set = Vector{Int64}(undef, r)
                entries_set[1:i] .= -1
                entries_set[(i + 1):(r - i + 1)] .= 0
                entries_set[(r - i + 2):r] .= 1
                
                # Add a generator of all permutations given `i` -1's to `generator_slices`
                generator_slices[idx] = (
                    vcat(leading, body) for body in multiset_permutations(entries_set, r)
                )
                idx += 1
            end
        end
        
        # Chain all generators together to form a single lazy iterator
        eigvec_generator = Iterators.flatten(generator_slices)
        return eigvec_generator
    end
end