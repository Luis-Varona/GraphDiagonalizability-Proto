module UniqueWithTolerance
    """
        uniquetol(vec; atol=1e-8, rtol=sqrt(ϵ),
                  return_indices=false, return_inverse=false,
                  return_counts=false, occurrence="highest")
    
    ADD LATER
    
    # Arguments
    - `vec::AbstractVector{T}`: ADD LATER
    
    # Keyword Arguments
    - `atol::Real`: ADD LATER
    - `rtol::Real`: ADD LATER
    - `return_indices::Bool`: ADD LATER
    - `return_inverse::Bool`: ADD LATER
    - `return_counts::Bool`: ADD LATER
    
    # Returns
    - `vec_unique::Vector{T}`: ADD LATER
    - `indices_unique::Vector{Int64}`: ADD LATER
    - `inverse_unique::Vector{Int64}`: ADD LATER
    - `counts_unique::Vector{Int64}`: ADD LATER
    
    # Examples
    ADD LATER
    """
    function uniquetol(
        vec::AbstractVector{T};
        atol::Real=1e-8,
        rtol::Real=sqrt(eps(real(float(one(T))))),
        return_indices::Bool=false,
        return_inverse::Bool=false,
        return_counts::Bool=false,
        occurrence::String="highest",
    ) where T
        if !(occurrence in ("highest", "lowest"))
            throw(ArgumentError("`occurrence` must be either \"highest\" or \"lowest\""))
        end
        
        n = length(vec)
        
        if n == 0
            vec_unique = T[]
            return_indices && (indices_unique = Int64[])
            return_inverse && (inverse_unique = Int64[])
            return_counts && (counts_unique = Int64[])
        else
            isclose(x, y) = isapprox(x, y; atol=atol, rtol=rtol)
            
            perm_sorted = sortperm(vec)
            sorted_vec = vec[perm_sorted]
            indices_unique = Int64[]
            idx = 1
            idx_switch = 0
            
            while !isnothing(idx_switch)
                idx += idx_switch
                push!(indices_unique, idx)
                c = sorted_vec[idx]
                idx_switch = findfirst(x -> !isclose(c, x), sorted_vec[(idx + 1):end])
            end

            k = length(indices_unique)
            
            if return_counts || return_inverse
                counts_unique = diff(vcat(indices_unique, n + 1))
                if return_inverse
                    inverse_unique = vcat(fill.(eachindex(counts_unique), counts_unique)...)
                    permute!(inverse_unique, invperm(perm_sorted))
                end
            end
            
            if occurrence == "highest"
                indices_unique[2:end] .-= 1
                indices_unique[1] = n
                permute!(indices_unique, vcat(2:k, 1))
            end
            
            indices_unique .= perm_sorted[indices_unique]
            vec_unique = Vector{T}(undef, k)
            vec_unique .= vec[indices_unique]
        end
        
        if return_indices || return_inverse || return_counts
            output_temp = Vector{Int64}[]
            return_indices && push!(output_temp, indices_unique)
            return_inverse && push!(output_temp, inverse_unique)
            return_counts && push!(output_temp, counts_unique)
            output = tuple(vec_unique, output_temp...)
        else
            output = vec_unique
        end
        
        return output
    end
    
    """
        uniquetol(vec, rtol;
                  return_indices=false, return_inverse=false,
                  return_counts=false, occurrence="highest")
    
    ADD LATER
    
    # Examples
    ADD LATER
    """
    uniquetol(
        vec::AbstractVector, rtol::Real;
        return_indices::Bool=false,
        return_inverse::Bool=false,
        return_counts::Bool=false,
        occurrence::String="highest",
    ) = uniquetol(
        vec;
        rtol=rtol,
        return_indices=return_indices,
        return_inverse=return_inverse,
        return_counts=return_counts,
        occurrence=occurrence,
    )
end