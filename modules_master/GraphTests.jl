module GraphTests
    #- EXPORTS AND IMPORTS
    using DataStructures: Stack
    using Graphs
    using MolecularGraph: nodesubgraph_is_isomorphic
    
    
    #- FUNCTION: `is_bipartite`
    """
        is_bipartite(G)
    
    Determine whether an undirected graph is bipartite.
    
    # Arguments
    - `G::SimpleGraph`: an undirected graph.
    
    # Keyword Arguments
    - `return_bipartition::Bool=false`: whether to return the bipartition of the graph.
    
    # Returns
    ADD LATER
    """
    function is_bipartite(A::AbstractMatrix; return_bipartition::Bool=false)
        return is_bipartite(SimpleGraph(A); return_bipartition=return_bipartition)
    end
    
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
        
        return return_bipartition ? (bipartite, coloring) : bipartite
    end
    
    
    #- FUNCTION: `is_cograph`
    """
        is_cograph(G)
    
    Determine whether an undirected graph is a cograph.
    
    # Arguments
    - `G::SimpleGraph`: an undirected graph.
    
    # Returns
    - `Bool`: whether `G` is a cograph.
    """
    function is_cograph(G::SimpleGraph)
        return !nodesubgraph_is_isomorphic(G, path_graph(4))
    end
end