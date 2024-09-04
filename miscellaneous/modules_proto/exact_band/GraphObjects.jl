module GraphObjects
    #- EXPORTS AND IMPORTS
    export DiagGraph, DiagGraphs_to_PyArray, graph6_to_laplacian
    
    using LinearAlgebra: Diagonal, diag, I
    using Intervals: Interval, Closed
    using PythonCall: Py, pyimport, pylist, pytuple, pyconvert
    np = pyimport("numpy")
    nx = pyimport("networkx")
    
    
    #- STRUCT: `DiagGraph`
    """
        DiagGraph{Int64} <: AbstractMatrix{Int64}
    
    An undirected graph in Laplacian form with data on its `{-1,0,1}`- and `{-1,1}`-spectra.
    
    # Fields
    - `laplacian_matrix::Matrix{Int64}`: the (integer-weighted) Laplacian matrix.
    - `band_zerooneneg::Union{Float64, Int64, Interval{Float64, Closed, Closed}}`: the
        `{-1,0,1}`-bandwidth of the graph. If only tested up to a certain
        `{-1,0,1}`-bandwidth `k` (e.g., `2` for weak Hadamard diagonalizability) and the
        result is `false`, this is set to a closed interval between `k + 1` and `Inf`.
    - `band_oneneg::Union{Float64, Int64, Interval{Float64, Closed, Closed}}`: the
        `{-1,1}`-bandwidth of the graph. If only tested up to a certain `{-1,1}`-bandwidth
        `k` (e.g., `1` for Hadamard diagonalizability) and the result is `false`, this is
        set to a closed interval between `k + 1` and `Inf`.
    - `eigvals::Vector{Int64}`: the eigenvalues of the Laplacian matrix. If the graph is not
        Laplacian integral (and thus not `{-1,0,1}`-diagonalizable), an empty array.
    - `eigvecs_zerooneneg::Matrix{Int64}`: a `{-1,0,1}`-matrix of Laplacian eigenvectors. If
        the graph is not `{-1,0,1}`-diagonalizable, an empty array.
    - `eigvecs_oneneg::Matrix{Int64}`: a `{-1,1}`-matrix of Laplacian eigenvectors. If the
        graph is not `{-1,1}`-diagonalizable, an empty array.
    
    # Methods
    - `order::Int64`: the number of vertices in the graph.
    - `size::Int64`: the number of edges in the graph.
    - `density::Float64`: the (unweighted) density of the graph.
    - `weighted_density::Float64`: the weighted density of the graph.
    - `average_degree::Float64`: the (unweighted) average degree of the graph.
    - `average_weighted_degree::Float64`: the weighted average degree of the graph.
    
    - `adjacency_matrix::Matrix{Int64}`: the adjacency matrix of the graph.
    - `graph6_string::String`: the graph6 string of the graph.
    
    - `is_weighted::Bool`: whether the graph has edge weights.
    - `is_negatively_weighted::Bool`: whether the graph has negative edge weights.
    - `is_regular::Bool`: whether the graph is regular.
    - `is_bipartite::Bool`: whether the graph is bipartite.
    - `is_planar::Bool`: whether the graph is planar.
    """
    struct DiagGraph
        laplacian_matrix::Matrix{Int64}
        band_zerooneneg::Union{Float64, Int64, Interval{Float64, Closed, Closed}}
        band_oneneg::Union{Float64, Int64, Interval{Float64, Closed, Closed}}
        eigvals::Vector{Int64}
        eigvecs_zerooneneg::Matrix{Int64}
        eigvecs_oneneg::Matrix{Int64}
    end
    

    #- String representation of `DiagGraph`
    function Base.show(io::IO, g::DiagGraph)
        if g.band_zerooneneg < Inf
            print(io, "{-1,0,1}-diagonalizable graph on $(g.order) vertices")
        else
            print(io, "{-1,0,1}-non-diagonalizable graph on $(g.order) vertices")
        end
    end
    
    
    #- Outer constructor for `DiagGraph`
    function Base.getproperty(g::DiagGraph, prop::Symbol)
        attributes = (
            :laplacian_matrix,
            :band_zerooneneg,
            :band_oneneg,
            :eigvals,
            :eigvecs_zerooneneg,
            :eigvecs_oneneg
        )
        
        if prop in attributes
            return getfield(g, prop)
        end
        
        if prop == :order
            return size(g.laplacian_matrix, 1)
        end
        if prop == :size
            A_unweighted = (g.adjacency_matrix .!= 0)
            return Int64(sum(A_unweighted) / 2)
        end
        if prop == :density
            n = g.order
            s = g.size
            return s == 1 ? NaN : 2s / (n * (n - 1))
        end
        if prop == :weighted_density
            n = g.order
            s = g.size
            return s == 1 ? NaN : sum(g.adjacency_matrix) / (n * (n - 1))
        end
        if prop == :average_degree
            return g.size / g.order
        end
        if prop == :average_weighted_degree
            return sum(g.adjacency_matrix) / (2g.order)
        end
        
        if prop == :adjacency_matrix
            L = g.laplacian_matrix
            return Diagonal(diag(L)) - L
        end
        if prop == :graph6_string
            return networkx_to_graph6(g.networkx_graph)
        end
        
        if prop == :is_weighted
            return !all(x -> x in [0, 1], g.adjacency_matrix)
        end
        if prop == :is_negatively_weighted
            return any(g.adjacency_matrix .< 0)
        end
        if prop == :is_regular
            return pyconvert(Bool, nx.is_regular(
                nx.from_numpy_array(np.array(g.adjacency_matrix)))
            )
        end
        if prop == :is_bipartite
            return pyconvert(Bool, nx.is_bipartite(
                nx.from_numpy_array(np.array(g.adjacency_matrix)))
            )
        end
        if prop == :is_planar
            return pyconvert(Bool, nx.is_planar(
                nx.from_numpy_array(np.array(g.adjacency_matrix)))
            )
        end
    end
    
    
    #- FUNCTION: 'DiagGraphs_to_PyArray'
    """
        DiagGraphs_to_PyArray(graphs::Vector{DiagGraph})::Py
    
    Convert a vector of `DiagGraph` instances to a NumPy object array.
    
    # Arguments
    - `graphs::Vector{DiagGraph}`: a vector of `DiagGraph` instances.
    
    # Returns
    - `Py`: a NumPy object array of `graphs`.
    
    # Examples
    ADD LATER
    """
    function DiagGraphs_to_PyArray(graphs::Vector{DiagGraph})::Py
        function DiagGraph_to_PyTuple(g::DiagGraph)::Py
            return pytuple([
                g.graph6_string,
                g.band_zerooneneg,
                g.band_oneneg,
                np.array(g.eigvals),
                np.array(g.eigvecs_zerooneneg),
                np.array(g.eigvecs_oneneg)
            ])
        end
        
        max_order = maximum([g.order for g in graphs])
        PyGraph = np.dtype(pylist([
            ("graph6_string", "<U$max_order"),
            ("band_zerooneneg", "<i8"),
            ("band_oneneg", "O"),
            ("eigvals", "O"),
            ("eigvecs_zerooneneg", "O"),
            ("eigvecs_oneneg", "O")
        ]))
        
        return np.array(pylist(DiagGraph_to_PyTuple.(graphs)), dtype = PyGraph)
    end
    

    #- FUNCTION: `graph6_to_laplacian`
    """
        graph6_to_laplacian(s::String)::Matrix{Int64}
    
    Convert a graph from graph6 format to a Laplacian matrix.
    
    # Arguments
    - `s::String`: a graph6 string.
    
    # Returns
    - `Matrix{Int64}`: the Laplacian matrix of the graph.
    
    # Examples
    Produce the Laplacian of the Cartesian product K_2 □ K_3:
    ```jldoctest
    julia> s = "E{Sw"
    julia> graph6_to_laplacian(s)
    6×6 Matrix{Int64}:
      3  -1  -1  -1   0   0
     -1   3  -1   0  -1   0
     -1  -1   3   0   0  -1
     -1   0   0   3  -1  -1
      0  -1   0  -1   3  -1
      0   0  -1  -1  -1   3
    """
    function graph6_to_laplacian(s::String)::Matrix{Int64}
        g = nx.from_graph6_bytes(Base.CodeUnits(s))
        A = pyconvert(Matrix{Int64}, nx.to_numpy_array(g))
        return I(size(A, 1)) .* sum(A, dims = 2) - A
    end
    
    
    #- FUNCTION: `networkx_to_graph6`
    """
        networkx_to_graph6(g::Py)::String
    
    Convert a NetworkX (Python) graph to a graph6 Julia string.
    
    # Arguments
    - `g::Py`: a NetworkX graph
    
    # Returns
    - `String`: the graph6 string of `g`
    """
    function networkx_to_graph6(g::Py)::String
        return String(pyconvert(Base.CodeUnits, nx.to_graph6_bytes(g)))[11:(end - 1)]
    end
end