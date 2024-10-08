module KOrthogonalizability
    #- EXPORTS AND IMPORTS
    export is_k_orthogonalizable    
    
    using LinearAlgebra: diagind
    using Graphs: SimpleGraph
    using MolecularGraph: subgraph_monomorphisms
    
    include("OrthogonalityGraphs.jl")
    using .OrthogonalityGraphs: orthograph_complement
    
    
    #- FUNCTION: `is_k_orthogonalizable`
    """
        is_k_orthogonalizable(X::AbstractMatrix, k::Int64)
    
    Determine whether a matrix is `k`-orthogonalizable and find a column ordering if so.
    
    Multiple dispatch on `orthograph_complement` is used to determine whether exact or
    approximate orthogonality between columns is required.
    
    # Arguments
    - `X::AbstractMatrix`: the matrix to test for `k`-orthogonalizability.
    - `k::Int64`: the orthogonalizability parameter.
    
    # Returns
    - `is_k_ortho`::Bool: whether `X` is k-orthogonalizable.
    - `order::Vector{Int64}`: if `is_k_ortho`, a `k`-orthogonal column ordering of `X`.
        Otherwise, an empty array.
    
    # Examples
    Test an integer matrix with minimum `k`-orthogonalizability `4`:
    ```jldoctest
    julia> X = [ -2   4   -1    0   12;
                  7  -1  -13   -8  -13;
                 -2   1    8    7    8;
                 -4   4   -5  -11   -9;
                 12  -8    1   10    3;
                -13  -1   -7    0   13;
                  9  -2   -2   12    4];
    julia> is_k_orthogonalizable(X, 1)
    (false, Int64[])
    ```
    
    Test a (complex) integer matrix with minimum `k`-orthogonalizability `3`:
    ```jldoctest
    julia> X = [-5-2im  0+0im    2+6im   -4+14im;
                -3-1im  0+0im   -1+11im   6+6im;
                -2+2im  4-10im   0+0im    1+15im];
    julia> is_k_orthogonalizable(X, 3)
    (true, [2, 4, 1, 3])
    ```
    
    Test a floating-point matrix with minimum `k`-orthogonalizability `2`:
    ```jldoctest
    julia> X = [-1.864046107   1.21905858    2.343251229   0.955248594;
                 1.118427664  -2.438117161  -7.029753687   0.955248594;
                -1.864046107  -8.533410062   4.686502458  -0.955248594;
                -2.236855328   4.876234321  -9.373004916  -1.432872891];
    julia> is_k_orthogonalizable(X, 3)
    (true, [1, 4, 2, 3])
    ```
    """
    function is_k_orthogonalizable(X::AbstractMatrix, k::Int64)::Tuple{Bool, Vector{Int64}}
        n = size(X, 2)
        
        if k >= n
            is_k_ortho = true # Every m × k matrix is k-orthogonal
            order = collect(1:n) # Any matrix, and hence any column ordering, works
        elseif k == 1 # Equivalent to pairwise orthogonality
            # Pairwise orthogonality is permutation-invariant
            is_k_ortho = !any(orthograph_complement(X))
            order = is_k_ortho ? collect(1:n) : Int64[]
        elseif k == 2 # Equivalent to quasi-orthogonalizability
            (is_k_ortho, order) = is_quasi_orthogonalizable(X)
        else
            H = orthograph_complement(X)
            
            if any(sum(H, dims = 2) .> 2k - 2)
                is_k_ortho = false # A column can be non-orthogonal to at most 2k - 2 others
                order = Int64[]
            else
                # Construct the adjacency matrix of the Turtle graph T(n, k)
                G = zeros(Bool, n, n)
                for i in vcat((1 - k):-1, 1:(k - 1)) # Make all possible columns non-orthogonal
                    G[diagind(G, i)] .= true
                end
                
                # Determine whether `X` has at least the minimum set of orthogonal pairs
                monomorphism = iterate(
                    subgraph_monomorphisms(SimpleGraph(G), SimpleGraph(H))
                )
                
                if isnothing(monomorphism)
                    is_k_ortho = false
                    order = Int64[]
                else
                    is_k_ortho = true
                    # The monomorphism mapping determines the k-orthogonal column ordering
                    order = last.(sort(collect(monomorphism[1])))
                    
                    if order[1] > order[end]
                        order = reverse(order) # To better preserve the original ordering
                    end
                end
            end
        end
        
        return (is_k_ortho, order)
    end
    
    
    #- FUNCTION: `is_quasi_orthogonalizable`
    """
        is_quasi_orthogonalizable(X::AbstractMatrix)
    
    Determine whether a matrix is quasi-orthogonalizable and find a column ordering if so.
    
    # Arguments
    - `X::AbstractMatrix`: the matrix to test for quasi-orthogonalizability.
    
    # Returns
    - `is_quasi::Bool`: whether `X` is quasi-orthogonalizable.
    - `order::Vector{Int64}`: if `is_quasi`, a quasi-orthogonal column ordering of `X`.
        Otherwise, an empty array.
    """
    function is_quasi_orthogonalizable(X::AbstractMatrix)::Tuple{Bool, Vector{Int64}}
        n = size(X, 2)
        
        if n <= 2
            is_quasi = true # Every m × 2 matrix is quasi-orthogonal
            order = collect(1:n) # Any matrix, and hence any column ordering, works
        else
            A = orthograph_complement(X)
            vertex_degrees = sum.(eachrow(A))
            
            if all(vertex_degrees .== 0)
                is_quasi = true # In fact, `X` is pairwise orthogonal
                order = collect(1:n) # Pairwise orthogonality is permutation-invariant
            elseif any(vertex_degrees .> 2)
                is_quasi = false # A column can only be non-orthogonal to contiguous columns
                order = Int64[]
            elseif sum(vertex_degrees .== 1) < 2
                is_quasi = false # The orthogonality graph complement is not a path subgraph
                order = Int64[]
            else
                visited = zeros(Bool, n)
                order = Vector{Int64}(undef, n)
                
                # Identify the columns of `X` that are orthogonal to all others
                isolated_vertices = findall(x -> x == 0, vertex_degrees)
                visited[isolated_vertices] .= true
                num_isolated = length(isolated_vertices)
                
                # These columns can come first in the k-orthogonal column ordering
                order[1:num_isolated] = isolated_vertices
                idx_order = num_isolated + 1
                
                """
                    traverse_path!(visited, order, idx_order, vertex)
                
                Traverse a path starting at some `vertex` with one outgoing edge.
                
                Every vertex that can be visited by traversing the path is marked in
                `visited`. The order of visited vertices is recorded in `order`. This
                function only works when `vertex` has exactly one outgoing edge and the
                maximum vertex degree of `A` is 2 (as ensured in the outer function).
                
                # Arguments
                - visited::Vector{Bool}: the labels of previously visited vertices.
                - order::Vector{Int64}: the order in which vertices are visited. This
                    determines the k-orthogonal column ordering of `X`.
                - idx_order::Int64: the number of visited vertices + 1.
                - vertex::Int64: the label of the initial path endpoint.
                
                # Returns
                - visited::Vector{Bool}: the mutated `visited` vector.
                - order::Vector{Int64}: the mutated `order` vector.
                - idx_order::Int64: the updated `idx_order` value.
                """
                function traverse_path!(
                    visited::Vector{Bool},
                    order::Vector{Int64},
                    idx_order::Int64,
                    vertex::Int64
                )::Tuple{Vector{Bool}, Vector{Int64}, Int64}
                    
                    visited[vertex] = true
                    order[idx_order] = vertex
                    idx_order += 1
                    reached_endpoint = false
                    
                    while !reached_endpoint # Continue until the other endpoint is reached
                        for neighbor in 1:n # Search for the next neighbor of `vertex`
                            if !visited[neighbor] && A[vertex, neighbor]
                                visited[neighbor] = true
                                order[idx_order] = neighbor
                                idx_order += 1
                                vertex = neighbor # Move to `neighbor` in the path traversal
                                break
                            end
                        end
                        
                        # If the next vertex has degree 1, the path traversal is complete
                        reached_endpoint = (vertex_degrees[vertex] == 1)
                    end
                    
                    return (visited, order, idx_order)
                end
                
                for vertex in 1:n # Traverse all paths originating from degree-1 vertices
                    if !visited[vertex] && vertex_degrees[vertex] == 1
                        (visited, order, idx_order) = traverse_path!(
                            visited, order, idx_order, vertex
                        )
                    end
                end
                
                if order[1] > order[end]
                    order = reverse(order) # To better preserve the original ordering
                end
                
                if !all(visited)
                    is_quasi = false # Some connected component of the graph of `A` is a cycle
                    order = Int64[]
                else
                    is_quasi = true # Every connected component is a path or an isolated vertex
                end
            end
        end
        
        return (is_quasi, order)
    end
end