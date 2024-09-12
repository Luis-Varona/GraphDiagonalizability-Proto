module UniqueWithTolerance    
    function uniquetol(
        A::AbstractVector;
        atol::Real=1e-8,
        rtol::Real=sqrt(eps(real(float(one(eltype(A)))))),
        return_index::Bool=false,
        return_counts::Bool=false,
        occurrence::String="highest",
    )
        if !(occurrence in ("highest", "lowest"))
            throw(ArgumentError("`occurrence` must be either \"highest\" or \"lowest\""))
        end
        
        splits = _splittol(A; atol=atol, rtol=rtol)
    end
    
    function _uniquetol_when_closechain(
        A::AbstractVector;
        atol::Real, rtol::Real, return_index::Bool, return_counts::Bool, occurrence::String,
    )
        n = length(A)
        idxs = Int64[]
        idx = 1
        foo = 0
        
        while !isnothing(foo)
            idx += foo
            push!(idxs, idx)
            c = A[idx]
            foo = findfirst(x -> !isapprox(c, x; atol=atol, rtol=rtol), A[(idx + 1):end])
        end
        
        return_counts && (counts = diff(vcat(idxs, n + 1)))
        
        if occurrence == "highest"
            idxs[2:end] .-= 1
            idxs[1] = n
            permute!(idxs, vcat(2:length(idxs), 1))
        end
        
        unique_arr = A[idxs]
        
        output = if return_index && return_counts
            (unique_arr, idxs, counts)
        elseif return_index
            (unique_arr, idxs)
        elseif return_counts
            (unique_arr, counts)
        else
            unique_arr
        end
        
        return output
    end
    
    function _splittol(A::AbstractVector{T}; atol::Real, rtol::Real) where T
        n = length(A)
        a = sort(A)
        split_lengths = Int64[]
        idx = 1
        
        while idx < n
            ct = 0
            close = true
            
            while close && (idx < n) # Use `findfirst` instead?
                close = isapprox(a[idx], a[idx + 1]; atol=atol, rtol=rtol)
                idx += 1
                ct += 1
            end
            
            push!(split_lengths, ct)
        end
        
        (idx == n) && push!(split_lengths, 1)
        k = length(split_lengths)
        splits = [Vector{T}(undef, l) for l in split_lengths]
        start = 1
        
        for i in 1:k
            stop = start + split_lengths[i] - 1
            splits[i] .= a[start:stop]
            start = stop + 1
        end
        
        return splits
    end
    
    
    function _closechain(A::AbstractVector; atol::Real, rtol::Real)
        foo(idx) = isapprox(A[idx], A[idx + 1]; atol=atol, rtol=rtol)
        return all(foo.(1:(length(A) - 1)))
    end
    
    # rtol::Real=1e-5
    # rtol::Real=iszero(atol) * length(A) * eps(real(float(one(eltype(A)))))
    
    
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
    
    
    a = [1 + x*1e-9 for x in -107:108]
end