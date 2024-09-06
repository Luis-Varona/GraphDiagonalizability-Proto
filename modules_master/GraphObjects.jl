module GraphObjects
    #- EXPORTS AND IMPORTS
    export DiagGraph
    using Intervals: Interval, Closed
    
    #- TYPE: `DiagGraph`
    struct DiagGraph{
        β<:Union{Float64, Int64, Interval{Int64, Closed, Closed}},
        B<:Union{Float64, Int64, Interval{Int64, Closed, Closed}},
        λ<:Union{Missing, Int64},
        τ<:Union{Missing, λ},
        T<:Union{Missing, τ},
    }
        laplacian_matrix::Matrix{Int64}
        band_zerooneneg::β
        band_oneneg::B
        eigvals::Vector{λ}
        eigvecs_zerooneneg::Matrix{τ}
        eigen_oneneg::Matrix{T}
    end
end