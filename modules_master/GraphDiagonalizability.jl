module GraphDiagonalizability
    #- EXPORTS AND IMPORTS
    export bandwidths_zerooneneg
    
    include("GraphObjects.jl")
    include("KOrthogonalizability.jl")
    include("LaplacianSpectra.jl")
    using .GraphObjects: DiagGraph
    using .KOrthogonalizability: is_k_orthogonalizable
    using .LaplacianSpectra: eigvecs_zerooneneg, is_spectrum_integral
    
    using Combinatorics: combinations
    using Intervals: Interval, Closed
    using LinearAlgebra: rank
    
    
    #- FUNCTION: `bandwidths_zerooneneg`
    """
        bandwidths_zerooneneg(
            L;
            min_zerooneneg=1,
            max_zerooneneg=size(L, 1),
            min_oneneg=1,
            max_oneneg=size(L, 1),
        )
    
    ADD LATER
    
    # Arguments
    - `L::AbstractMatrix{Int64}`: ADD LATER

    # Keyword Arguments
    - `min_zerooneneg::Int64=1`: ADD LATER
    - `max_zerooneneg::Int64=size(L, 1)`: ADD LATER
    - `min_oneneg::Int64=1`: ADD LATER
    - `max_oneneg::Union{Float64, Int64}=size(L, 1)`: ADD LATER
    
    # Returns
    ADD LATER
    
    # Examples
    ADD LATER
    """
    function bandwidths_zerooneneg(
        L::AbstractMatrix{Int64};
        min_zerooneneg::Int64=1,
        max_zerooneneg::Int64=size(L, 1),
        min_oneneg::Union{Float64, Int64}=1,
        max_oneneg::Int64=size(L, 1),
    )
        n = size(L, 1)
        (spectrum_integral, λs_sorted, λ_counts) = is_spectrum_integral(L)
        
        if !spectrum_integral
            Γ = DiagGraph(L, Inf, Inf, missing, missing, missing)
        else
            (
                eigvecs,
                cols_oneneg,
                bases_zerooneneg_temp,
                bases_oneneg_temp,
            ) = eigvecs_zerooneneg(L, λ_counts)
            
            multis = last.(λ_counts)
            
            if any(size.(bases_zerooneneg_temp, 2) .<= multis)
                Γ = DiagGraph(L, Inf, Inf, λs_sorted, missing, missing)
            else
                r = length(λ_counts)
                sup_zerooneneg = max(multis)
                initial_bands = fill(sup_zerooneneg, r)
                
                for k in (sup_zerooneneg - 1):-1:1
                    results = is_k_orthogonalizable.(bases_zerooneneg, k)
                    initial_bands[results] = k
                    
                    all(results) && (sup_zerooneneg = k; break)
                    any(results) || break
                end
                
                if any(size.(bases_oneneg, 2) .<= multis)
                    min_oneneg = Inf
                else
                    sup_oneneg = sup_zerooneneg
                end
                
                bases_zerooneneg = Vector{Matrix{Int64}}(undef, r)
                bases_oneneg = Vector{Matrix{Int64}}(undef, r)
                
                for ((idx, band), (λ, μ)) in zip(enumerate(initial_bands), λ_counts)
                    (
                        min_zerooneneg,
                        min_oneneg,
                        bases_zerooneneg[idx],
                        bases_oneneg[idx],
                    ) = _eigspace_bandwidths(
                        L, λ, μ, eigvecs[idx], cols_oneneg[idx];
                        min_zerooneneg=min_zerooneneg,
                        max_zerooneneg=min(band - 1, max_zerooneneg),
                        min_oneneg=min_oneneg,
                        max_oneneg=max_oneneg,
                    )
                    
                    if isa(min_zerooneneg, Interval)
                        min_zerooneneg = Interval{Closed, Closed}(
                            min_zerooneneg.first,
                            sup_zerooneneg - 1,
                        )
                        break
                    end
                    
                    (min_zerooneneg > max_zerooneneg) && (P_zerooneneg = bases_zerooneneg_temp)
                    (Inf > min_oneneg > max_oneneg) && (P_oneneg = bases_oneneg_temp)
                    (min_zerooneneg > max_zerooneneg) && (min_oneneg > max_oneneg) && break
                end
                
                (min_zerooneneg <= max_zerooneneg) && (P_zerooneneg = hcat(bases_zerooneneg...)) # ADD LATER (type instability)
                (min_oneneg <= max_oneneg) && (P_oneneg = hcat(bases_oneneg...)) # ADD LATER (type instability)
                (min_oneneg == Inf) && (P_oneneg = missing) # ADD LATER (type instability)
                
                Γ = DiagGraph(L, min_zerooneneg, min_oneneg, λs_sorted, P_zerooneneg, P_oneneg)
            end
        end
        
        return Γ
    end
    
    
    #- FUNCTION: `_eigspace_bandwidths`
    """
        ADD LATER
    """
    function _eigspace_bandwidths(
        L::AbstractMatrix{Int64},
        λ::Int64,
        μ::Int64,
        eigspace::Matrix{Int64},
        cols_oneneg::Vector{Int64};
        min_zerooneneg::Int64,
        max_zerooneneg::Int64,
        min_oneneg::Union{Float64, Int64},
        max_oneneg::Int64,
    )
        n = size(L, 1)
        
        if μ == 1
            band_zerooneneg = min_zerooneneg
            band_oneneg = min_oneneg
            
            if band_oneneg < Inf
                basis_zerooneneg = eigspace[:, cols_oneneg[1]]
                basis_oneneg = copy(basis_zerooneneg)
            else
                basis_zerooneneg = eigspace[:, 1]
                basis_oneneg = Matrix{Int64}(undef, n, 0)
            end
        else
            function DFS(idxs::Vector{Int64}, k::Int64, idx_set::Vector{Int64})
                depth = length(idxs)
                partial_basis = eigspace[:, idxs]
                
                if rank(partial_basis, 1e-5) < depth
                    k_ortho = false
                    idxs_final = Int64[]
                else
                    (k_ortho, order) = is_k_orthogonalizable(partial_basis, k)
                    
                    if !k_ortho
                        idxs_final = Int64[]
                    elseif depth == μ
                        idxs_final = idxs[order]
                    else
                        k_ortho = false
                        idxs_final = Int64[]
                        last = idxs[end]
                        
                        for idx_next in filter(i -> i > last, idx_set)
                            (k_ortho, idxs_final) = DFS(vcat(idxs, idx_next), k, idx_set)
                            k_ortho && break
                        end
                    end
                end
                
                return (k_ortho, idxs_final)
            end
            
            function search_roots(idx_set::Vector{Int64}, min_band::Int64, max_band::Int64)
                has_basis = false
                idxs_basis = Int64[]
                band = min_band
                
                while !has_basis && (band <= max_band)
                    for root in combinations(idx_set, band)
                        if rank(eigspace[:, root], 1e-5) == band
                            (has_basis, idxs_basis) = DFS(root, band, idx_set)
                            has_basis && break
                        end
                    end
                    
                    band += 1
                end
                
                basis = has_basis ? eigspace[:, idxs_basis] : Matrix{Int64}(undef, n, 0)
                return (band, basis)
            end
            
            (band_oneneg, basis_oneneg) = if min_oneneg == Inf
                (Inf, missing)
            else
                search_roots(cols_oneneg, min_oneneg, max_oneneg)
            end
            
            (band_zerooneneg, basis_zerooneneg) = if band_oneneg <= min_zerooneneg
                (min_zerooneneg, copy(basis_oneneg))
            else
                search_roots(
                    1:size(eigspace, 2),
                    min_zerooneneg,
                    min(band_oneneg - 1, max_zerooneneg),
                )
            end
            (band_zerooneneg >= band_oneneg) && (basis_zerooneneg = copy(basis_oneneg))
        end
        
        return (band_zerooneneg, band_oneneg, basis_zerooneneg, basis_oneneg)
    end
end