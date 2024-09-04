module WeakHadamards
    export get_weak_hadamards
    
    using LinearAlgebra: rank
    using Combinatorics: combinations
    
    include("../modules/KOrthogonalizability.jl")
    using .KOrthogonalizability: is_k_orthogonalizable
    
    function get_weak_hadamards(n::Int64)
        num_eigvecs = convert(Int64, (3^n - 1)/2)
        V = Vector{Vector{Int64}}(undef, num_eigvecs)
        i = 1
        
        for k in 1:n
            vec_slice1 = vcat(zeros(Int64, k - 1), 1)
            eigvecs_iter = Iterators.product(
                Iterators.repeated((-1, 0, 1), n - k)...)
            for vec_slice2 in eigvecs_iter
                V[i] = vcat(vec_slice1, collect(vec_slice2))
                i += 1
            end
        end
        
        if n <= 2
            idxs_bases = combinations(1:num_eigvecs, n)
        else
            idxs_start_opts = combinations(1:num_eigvecs, 2)
            depth = 3
            
            while depth <= n
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
        
        bases = [hcat(V[idxs]...) for idxs in idxs_bases]
        return bases
    end
end