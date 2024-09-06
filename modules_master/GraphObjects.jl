module GraphObjects
    #- EXPORTS AND IMPORTS
    export DiagGraph
    
    using Intervals: Interval, Closed
    using LinearAlgebra: Eigen
    
    
    #- TYPE: `DiagGraph`
    struct DiagGraph
        laplacian_matrix::Matrix{Int64}
        band_zerooneneg::Union{Int64, Interval{Int64, Closed, Closed}}
        band_oneneg::Union{Int64, Interval{Int64, Closed, Closed}}
        eigen_zerooneneg::Eigen{T, T, Matrix{T}, Vector{T}} where T<:Union{Missing, Int64}
        eigen_oneneg::Eigen{τ, τ, Matrix{τ}, Vector{τ}} where τ<:Union{Missing, Int64}
    end
end