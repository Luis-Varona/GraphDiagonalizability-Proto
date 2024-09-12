module UniqueWithTolerance    
    uniquetol(
        vec::AbstractVector, rtol::Real;
        return_indices::Bool=false, return_counts::Bool=false, occurrence::String="highest",
    ) = uniquetol(
        vec;
        rtol=rtol,
        return_indices=return_indices,
        return_counts=return_counts,
        occurrence=occurrence,
    )
    
    function uniquetol(
        vec::AbstractVector{T};
        atol::Real=1e-8,
        rtol::Real=sqrt(eps(real(float(one(T))))),
        return_indices::Bool=false,
        return_counts::Bool=false,
        occurrence::String="highest",
    ) where T
        if !(occurrence in ("highest", "lowest"))
            throw(ArgumentError("`occurrence` must be either \"highest\" or \"lowest\""))
        end
        
        n = length(vec)
        
        if n == 0
            vals_unique = T[]
            idxs_unique = Int64[]
            counts_unique = Int64[]
        else
            isclose(x, y) = isapprox(x, y; atol=atol, rtol=rtol)
            
            perm_sorted = sortperm(vec)
            sorted_vec = permute!(vec, perm_sorted)
            idxs_unique = Int64[]
            idx = 1
            idx_switch = 0
            
            while !isnothing(idx_switch)
                idx += idx_switch
                push!(idxs_unique, idx)
                c = sorted_vec[idx]
                idx_switch = findfirst(x -> !isclose(c, x), sorted_vec[(idx + 1):end])
            end
            
            return_counts && (counts_unique = diff(vcat(idxs_unique, n + 1)))
            
            if occurrence == "highest"
                idxs_unique[2:end] .-= 1
                idxs_unique[1] = n
                permute!(idxs_unique, vcat(2:length(idxs_unique), 1))
            end
            
            idxs_unique .= perm_sorted[idxs_unique]
            vals_unique = vec[idxs_unique]
        end
        
        output = if return_indices && return_counts
            (vals_unique, idxs_unique, counts_unique)
        elseif return_indices
            (vals_unique, idxs_unique)
        elseif return_counts
            (vals_unique, counts_unique)
        else
            vals_unique
        end
        
        return output
    end
    
    # function bar(x, n)
    #     R = rand(n)
    #     for r in R
    #         x ^= r
    #     end
    #     for r in R
    #         x ^= (1/r)
    #     end
    #     return x
    # end
    
    # t = [3.0, 3.0, 3.0, 4.0, 4.0, 7.0]
    # N = [5, 5, 3, 7, 4, 6]
    # t = bar.(t, N)
    # a = [1 + x*1e-9 for x in -107:108]
end