module GraphObjects
    #- EXPORTS AND IMPORTS
    export DiagGraph
    
    include("GraphFunctions.jl")
    using .GraphFunctions: is_bipartite, is_cartesian_product, is_cograph
    
    using Intervals: Interval, Closed
    using LinearAlgebra: Diagonal, diag
    # using PythonCall: Py, pyconvert, pyimport, pylist, pytuple
    
    # np = pyimport("numpy")
    # nx = pyimport("networkx")
    
    
    #- TYPE: `_BandwidthInterval`
    _BandwidthInterval = Interval{Int64, Closed, Closed}
    
    
    #- TYPE: `DiagGraph`
    struct DiagGraph{
        β,
        B,
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
        
        function DiagGraph(
            laplacian_matrix::Matrix{Int64},
            band_zerooneneg::β,
            band_oneneg::B,
            eigvals::Vector{λ},
            eigvecs_zerooneneg::Matrix{τ},
            eigen_oneneg::Matrix{T}
        ) where {β<:Float64, B<:Float64, λ, τ, T}
            return new{β, B, λ, τ, T}(
                laplacian_matrix,
                band_zerooneneg,
                band_oneneg,
                eigvals,
                eigvecs_zerooneneg,
                eigen_oneneg,
            )
        end
        
        function DiagGraph(
            laplacian_matrix::Matrix{Int64},
            band_zerooneneg::β,
            band_oneneg::B,
            eigvals::Vector{λ},
            eigvecs_zerooneneg::Matrix{τ},
            eigen_oneneg::Matrix{T}
        ) where {β<:_BandwidthInterval, B<:Union{Float64, _BandwidthInterval}, λ, τ, T}
            return new{β, B, λ, τ, T}(
                laplacian_matrix,
                band_zerooneneg,
                band_oneneg,
                eigvals,
                eigvecs_zerooneneg,
                eigen_oneneg,
            )
        end
        
        function DiagGraph(
            laplacian_matrix::Matrix{Int64},
            band_zerooneneg::β,
            band_oneneg::B,
            eigvals::Vector{λ},
            eigvecs_zerooneneg::Matrix{τ},
            eigen_oneneg::Matrix{T}
        ) where {β<:Int64, B<:Union{Float64, Int64, _BandwidthInterval}, λ, τ, T}
            return new{β, B, λ, τ, T}(
                laplacian_matrix,
                band_zerooneneg,
                band_oneneg,
                eigvals,
                eigvecs_zerooneneg,
                eigen_oneneg,
            )
        end
    end
    
    
    #- String representation of `DiagGraph`
    function Base.show(io::IO, g::DiagGraph)
        diag_signal = (g.band_zerooneneg != Inf) ? "" : "-non"
        print(io, "{-1,0,1}$diag_signal-diagonalizable graph on $(g.order) vertices")
    end
    
    
    #- Methods on `DiagGraph`
    function Base.getproperty(g::DiagGraph, prop::Symbol)
        if prop in (
            :laplacian_matrix,
            :band_zerooneneg,
            :band_oneneg,
            :eigvals,
            :eigvecs_zerooneneg,
            :eigvecs_oneneg
        )
            output = getfield(g, prop)
        
        elseif prop == :adjacency_matrix
            L = g.laplacian_matrix
            output = Diagonal(L) - L
        # elseif prop == :graph6_string
        #     output = ADD LATER
        
        elseif prop == :order
            output = size(g.laplacian_matrix, 1)
        elseif prop == :size
            output = Int64(sum(g.adjacency_matrix .!= 0) / 2)
        elseif prop == :density
            n = g.order
            s = g.size
            output = (s == 1) ? NaN : 2s / (n * (n - 1))
        elseif prop == :weighted_density
            n = g.order
            output = (g.size == 1) ? NaN : sum(g.adjacency_matrix) / (n * (n - 1))
        elseif prop == :average_degree
            output = g.size / g.order
        elseif prop == :average_weighted_degree
            output = sum(g.adjacency_matrix) / (2g.order)
        
        elseif prop == :is_connected
            output = is_connected(SimpleGraph(g.adjacency_matrix))
        elseif prop == :is_weighted
            output = any(x -> !(x in [0, 1]), g.adjacency_matrix)
        elseif prop == :is_negatively_weighted
            output = any(g.adjacency_matrix .< 0)
        elseif prop == :is_regular
            output = allequal(sum(g.adjacency_matrix .!= 0, dims = 2))
        elseif prop == :is_weighted_regular
            output = allequal(diag(g.laplacian_matrix))
        
        elseif prop == :is_cograph
            output = is_cograph(g.adjacency_matrix)
        elseif prop == :is_bipartite
            output = is_bipartite(g.adjacency_matrix)
        elseif prop == :is_cart_product
            output = is_cartesian_product(g.adjacency_matrix)
        elseif prop == :is_cart_product_complement
            A_complement = .!g.adjacency_matrix
            A_complement[diagind(A_complement)] .= false
            output = is_cartesian_product(A_complement)
        end
        
        return output
    end
    
    
    # #- FUNCTION: `networkx_to_graph6`
    # function networkx_to_graph6(g::Py)
    #     return String(pyconvert(Base.CodeUnits, nx.to_graph6_bytes(g)))[11:(end - 1)]
    # end
end