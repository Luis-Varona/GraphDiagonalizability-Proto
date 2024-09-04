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
                zerooneneg_bases,
                oneneg_bases,
            ) = eigvecs_zerooneneg(L, λ_counts)
            
            if !all(zerooneneg_bases)
                Γ = DiagGraph(L, Inf, Inf, λs_sorted, missing, missing)
            else
                sup_zerooneneg = max(last.(λ_counts))
                !all(oneneg_bases) ? (min_oneneg = Inf) : (sup_oneneg = sup_zerooneneg)
                
                for (idx, (λ, μ)) in enumerate(λ_counts)
                    (
                        min_zerooneneg,
                        min_oneneg,
                        basis_zerooneneg,
                        basis_oneneg,
                    ) = _eigspace_bandwidths(
                        L, λ, μ, eigvecs[idx], cols_oneneg[idx];
                        min_zerooneneg=min_zerooneneg,
                        max_zerooneneg=max_zerooneneg,
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
                    if min_zerooneneg > max_zerooneneg
                        # Get sup_zerooneneg basis for everything
                    end
                end
            end
        end
    end
    
    
    #- FUNCTION: `_eigspace_bandwidths`
    """
        ADD LATER
    """
    function _eigspace_bandwidths(
        L::AbstractMatrix{Int64},
        λ::Int64,
        μ::Int64,
        eigvecs::Matrix{Int64},
        cols_oneneg::Vector{Int64};
        min_zerooneneg::Int64,
        max_zerooneneg::Int64,
        min_oneneg::Union{Float64, Int64},
        max_oneneg::Int64,
    )
        3
    end
end