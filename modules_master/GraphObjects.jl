module GraphObjects
    #- EXPORTS AND IMPORTS
    export DiagGraph
    using Intervals: Interval, Closed
    
    #- TYPE: `DiagGraph`
    struct DiagGraph{
        β<:Union{Float64, Int64, Interval{Int64, Closed, Closed}},
        B<:Union{Float64, Int64, Interval{Int64, Closed, Closed}},
        τ<:Union{Missing, Int64},
        T<:Union{Missing, τ},
        S<:Union{Missing, T},
    }
        laplacian_matrix::Matrix{Int64}
        band_zerooneneg::β
        band_oneneg::B
        eigvals::Vector{τ}
        eigvecs_zerooneneg::Matrix{T}
        eigen_oneneg::Matrix{S}
    end
end