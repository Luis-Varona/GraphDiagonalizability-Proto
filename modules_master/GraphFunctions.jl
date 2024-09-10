module GraphFunctions
    #- EXPORTS AND IMPORTS
    using DataStructures: Stack
    using Graphs
    using MolecularGraph: nodesubgraph_is_isomorphic
    
    
    #- FUNCTION: `is_bipartite`
    """
        is_bipartite(G; return_bipartition=false)
        is_bipartite(A; return_bipartition=false)
    
    Determine whether an undirected graph is bipartite.
    
    Both SimpleGraphs and adjacency matrices are supported input formats.
    
    # Arguments
    - `G::SimpleGraph`: an undirected graph.
    - `A::AbstractMatrix`: an adjacency matrix (if `G` is not provided).
    
    # Keyword Arguments
    - `return_bipartition::Bool=false`: whether to return the bipartition of the graph.
    
    # Returns
    - `bipartite::Bool`: whether the graph is bipartite.
    - `coloring::Tuple{Vector{Int64}, Vector{Int64}}`: if `bipartite`, a partitioning of the
        graph vertices into two vectors. Otherwise, a tuple of empty vectors. (Only returned
        if `return_bipartition=true`.)
    """
    function is_bipartite(G::SimpleGraph; return_bipartition::Bool=false)
        bipartite = true
        coloring = Dict{Int64, Bool}()
        node = 1
        
        while bipartite && (node <= nv(G))
            if degree(G, node) == 0
                coloring[node] = false
            elseif !haskey(coloring, node)
                coloring[node] = true
                queue = Stack{Int64}()
                push!(queue, node)
                
                while !isempty(queue)
                    u = pop!(queue)
                    
                    for v in neighbors(G, u)
                        if !haskey(coloring, v)
                            coloring[v] = !coloring[u]
                            push!(queue, v)
                        elseif coloring[u] == coloring[v]
                            bipartite = false
                            empty!(coloring)
                            empty!(queue)
                            break
                        end
                    end
                end
            end
            
            node += 1
        end
        
        if return_bipartition
            bipartition = (
                collect(keys(filter(!last, coloring))),
                collect(keys(filter(last, coloring))),
            )
        end
        
        return return_bipartition ? (bipartite, bipartition) : bipartite
    end
    
    function is_bipartite(A::AbstractMatrix; return_bipartition::Bool=false)
        return is_bipartite(SimpleGraph(A); return_bipartition=return_bipartition)
    end
    
    
    #- FUNCTION: `is_cartesian_product`
    """
        is_cartesian_product(G)
        is_cartesian_product(A)
    
    Determine whether an undirected graph is a Cartesian product.
    
    Both SimpleGraphs and adjacency matrices are supported input formats.
    
    # Arguments
    - `G::SimpleGraph`: an undirected graph.
    - `A::AbstractMatrix`: an adjacency matrix (if `G` is not provided).
    
    # Returns
    - `Bool`: whether the graph is a Cartesian product.
    """
    function is_cartesian_product(G::SimpleGraph)
        3 # ADD LATER
    end
    
    # function is_cartesian_product(A::AbstractMatrix)
    #     return is_cartesian_product(SimpleGraph(A))
    # end
    
    
    #- FUNCTION: `is_cograph`
    """
        is_cograph(G)
        is_cograph(A)
    
    Determine whether an undirected graph is a cograph.
    
    Both SimpleGraphs and adjacency matrices are supported input formats.
    
    # Arguments
    - `G::SimpleGraph`: an undirected graph.
    - `A::AbstractMatrix`: an adjacency matrix (if `G` is not provided).
    
    # Returns
    - `Bool`: whether the graph is a cograph.
    """
    function is_cograph(G::SimpleGraph)
        return !nodesubgraph_is_isomorphic(G, path_graph(4))
    end
    
    function is_cograph(A::AbstractMatrix)
        return is_cograph(SimpleGraph(A))
    end
    
    
    #- FUNCTION: `graph6_string`
    function graph6_string(G::SimpleGraph)
        n = nv(G)
        lead = Char(nv(G) + 63)
        
        # ADD LATER
    end
    
    
    #- FUNCTION: `_binary_rep`
    function _binary_rep(G::SimpleGraph)
        n = nv(G)
        total_length = Int64(n * (n - 1) / 2)
        
        # ADD LATER
    end
end