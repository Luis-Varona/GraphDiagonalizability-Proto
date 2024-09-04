module KOrthogonalizability
    #- EXPORTS AND IMPORTS
    export is_k_orthogonalizable    
    
    include("OrthogonalityGraphs.jl")
    using .OrthogonalityGraphs: orthograph_complement
    
    using Graphs: SimpleGraph
    using LinearAlgebra: diagind
    using MolecularGraph: subgraph_monomorphisms
    
    
    #- FUNCTION: `is_k_orthogonalizable`
    """
        is_k_orthogonalizable(X, k)
    
    Determine whether a matrix is `k`-orthogonalizable and find a column ordering if so.
    
    Multiple dispatch on `orthograph_complement` is used to determine whether exact or
    approximate orthogonality between columns is required.
    
    # Arguments
    - `X::AbstractMatrix`: the matrix to test for `k`-orthogonalizability.
    - `k::Int64`: the orthogonalizability parameter.
    
    # Returns
    - `k_ortho`::Bool: whether `X` is k-orthogonalizable.
    - `order::Vector{Int64}`: if `k_ortho`, a `k`-orthogonal column ordering of `X`.
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
    function is_k_orthogonalizable(X::AbstractMatrix, k::Int64)
        n = size(X, 2)
        
        if k >= n # Every m × k matrix is k-orthogonal
            k_ortho = true
            order = collect(1:n) # Any matrix, and hence any column ordering, works
        elseif k == 1
            k_ortho = !any(orthograph_complement(X))
            # Pairwise orthogonality is permutation-invariant
            order = k_ortho ? collect(1:n) : Int64[]
        elseif k == 2
            (k_ortho, order) = _is_quasi_orthogonalizable(X)
        else
            H = orthograph_complement(X)
            
            # Each column can only be non-orthogonal to at most 2k - 2 others
            if any(sum(H, dims = 2) .> 2k - 2)
                k_ortho = false
                order = Int64[]
            else
                # Construct the adjacency matrix of the Turtle graph T(n, k)
                G = falses(n, n)
                for i in vcat((1 - k):-1, 1:(k - 1)) # Possible non-orthogonal column pairs
                    G[diagind(G, i)] .= true
                end
                
                # Determine whether `X` has at least the minimum set of orthogonal pairs
                monomorphism = iterate(
                    subgraph_monomorphisms(SimpleGraph(G), SimpleGraph(H))
                )
                if isnothing(monomorphism) # `X` cannot be rearranged to be k-orthogonal
                    k_ortho = false
                    order = Int64[]
                else
                    k_ortho = true
                    # The monomorphism mapping determines the k-orthogonal column ordering
                    order = last.(sort(collect(monomorphism[1])))
                    # To better preserve the original column ordering of `X`
                    (order[1] > order[end]) && (order = reverse(order))
                end
            end
        end
        
        return (k_ortho, order)
    end
    
    
    #- FUNCTION: `_is_quasi_orthogonalizable`
    """
        _is_quasi_orthogonalizable(X)
    
    Determine whether a matrix is quasi-orthogonalizable and find a column ordering if so.
    
    # Arguments
    - `X::AbstractMatrix`: the matrix to test for quasi-orthogonalizability.
    
    # Returns
    - `quasi_ortho::Bool`: whether `X` is quasi-orthogonalizable.
    - `order::Vector{Int64}`: if `quasi_ortho`, a quasi-orthogonal column ordering of `X`.
        Otherwise, an empty array.
    """
    function _is_quasi_orthogonalizable(X::AbstractMatrix)
        A = orthograph_complement(X)
        vertex_degrees = sum.(eachrow(A))
        
        if all(vertex_degrees .== 0) # `X` is pairwise orthogonal
            quasi_ortho = true
            order = collect(1:size(X, 2)) # Pairwise orthogonality is permutation-invariant
        # `A` cannot represent a subgraph of a path in either of these cases
        elseif any(vertex_degrees .> 2) || (sum(vertex_degrees .== 1) < 2)
            quasi_ortho = false
            order = Int64[]
        else
            # `X` is quasi-orthogonalizable if and only if `A` represents a path subgraph
            (quasi_ortho, order) = _is_path_subgraph(A, vertex_degrees)
        end
        
        return (quasi_ortho, order)
    end
    
    
    #- FUNCTION: `_is_path_subgraph`
    """
        _is_path_subgraph(A, vertex_degrees)
    
    ADD LATER
    
    [Note: works only because `_is_quasi_orthogonalizable` validates maximum degree `2`]
    
    # Arguments
    - `A::BitMatrix`: ADD LATER
    - `vertex_degrees::Vector{Int64}`: ADD LATER
    
    # Returns
    - `path_subgraph::Bool`: ADD LATER
    - `order::Vector{Int64}`: ADD LATER
    """
    function _is_path_subgraph(A::BitMatrix, vertex_degrees::Vector{Int64})
        n = size(A, 1)
        visited = falses(n)
        order = Vector{Int64}(undef, n)
        
        # Identify the columns of `X` that are orthogonal to all others
        isolated_vertices = findall(iszero, vertex_degrees)
        visited[isolated_vertices] .= true
        num_isolated = length(isolated_vertices)
        
        # Such columns can come first in the k-orthogonal column ordering
        order[1:num_isolated] = isolated_vertices
        idx_order = num_isolated + 1
        
        """
            _traverse_path!(visited, order, idx_order, node)
        
        Traverse a path starting at some `node` with one outgoing edge.
        
        Every node that can be visited by traversing the path is marked in `visited`. The
        order in which nodes are visited is recorded in `order`.
        
        # Arguments
        - visited::BitVector: the labels of previously visited node.
        - order::Vector{Int64}: the order in which nodes are visited. This determines the
            quasi-orthogonal column ordering of `X`.
        - idx_order::Int64: the number of visited nodes plus `1`.
        - node::Int64: the label of the initial path endpoint.
        
        # Returns
        - visited::BitVector: the mutated `visited` vector.
        - order::Vector{Int64}: the mutated `order` vector.
        - idx_order::Int64: the updated `idx_order` value.
        """
        function _traverse_path!(
            visited::BitVector,
            order::Vector{Int64},
            idx_order::Int64,
            node::Int64
        )
            visited[node] = true
            order[idx_order] = node
            idx_order += 1
            
            reached_endpoint = false
            while !reached_endpoint # Continue until the other endpoint is reached
                for neighbor in 1:n # Search for the next neighbor of `node`
                    if !visited[neighbor] && A[node, neighbor]
                        visited[neighbor] = true
                        order[idx_order] = neighbor
                        idx_order += 1
                        node = neighbor # Move to `neighbor` in the path traversal
                        break
                    end
                end
                # If the next node has degree 1, the path traversal is complete
                reached_endpoint = (vertex_degrees[node] == 1)
            end
            
            return (visited, order, idx_order)
        end
        
        for node in 1:n # Traverse all paths originating from degree-1 nodes
            if !visited[node] && (vertex_degrees[node] == 1)
                idx_order = _traverse_path!(visited, order, idx_order, node)[3]
            end
        end
        
        if !all(visited) # At least one component of the graph represented by `A` is a cycle
            path_subgraph = false
            order = Int64[]
        else # Every connected component of `A` is either a path or an isolated vertex
            path_subgraph = true
            # To better preserve the original column ordering of `X`
            (order[1] > order[end]) && (order = reverse(order))
        end
        
        return (path_subgraph, order)
    end
end