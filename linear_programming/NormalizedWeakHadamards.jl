module NormalizedWeakHadamards
    export normalized_weak_hadamards
    
    using LinearAlgebra: rank
    using Combinatorics: combinations, multiset_permutations
    
    include("../modules/KOrthogonalizability.jl")
    using .KOrthogonalizability: is_k_orthogonalizable
    
    function normalized_WH_columns(n::Int64)
        # V = [[1 for _ in 1:n]]
        V = Vector{Int64}[] # Note 1: Should pre-compute number of elements
        
        for k in 1:(n - 1)
            vec_slice1 = vcat(zeros(Int64, k - 1), 1)
            r = n - k
            for i in 1:ceil(Int64, r/2)
                entries_set = vcat(
                    [-1 for _ in 1:i],
                    [0 for _ in 1:(r - 2*i + 1)],
                    [1 for _ in 1:(i - 1)])
                eigvecs_iter = multiset_permutations(entries_set, r)
                for vec_slice2 in eigvecs_iter
                    push!(V, vcat(vec_slice1, collect(vec_slice2)))
                    # Fix later: See Note 1 (for efficient memory allocation)
                end
            end
        end
        return V
    end
    
    function normalized_weak_hadamards(n::Int64)
        V = normalized_WH_columns(n)
        num_eigvecs = length(V)
        
        if n <= 3
            idxs_bases = combinations(1:num_eigvecs, (n - 1))
        else
            idxs_start_opts = combinations(1:num_eigvecs, 2)
            depth = 3
            
            while depth <= (n - 1)
                idxs_new_opts = Vector{Int64}[]
                
                for idxs_start in idxs_start_opts
                    last = idxs_start[depth - 1]
                    for idx_next in (last + 1):num_eigvecs
                        idxs_comb = vcat(idxs_start, idx_next)
                        eigbasis = hcat(V[idxs_comb]...)
                        if rank(eigbasis, 1e-5) != depth
                            continue
                        end
                        
                        is_quasi, order = is_k_orthogonalizable(eigbasis, 2)
                        if is_quasi
                            push!(idxs_new_opts, idxs_comb[order])
                        end
                    end
                end
                
                idxs_start_opts = idxs_new_opts
                depth += 1
            end
            
            idxs_bases = idxs_start_opts
        end
        
        bases = [hcat([1 for _ in 1:n], V[idxs]...) for idxs in idxs_bases]
        return bases
    end
end