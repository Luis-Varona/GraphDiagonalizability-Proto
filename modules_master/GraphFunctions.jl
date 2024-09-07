module GraphFunctions
    #- EXPORTS AND IMPORTS
    using DataStructures: Queue
    using Graphs: SimpleGraph
    using MolecularGraph: nodesubgraph_is_isomorphic
    
    
    #- FUNCTION: `is_bipartite`
    """
        is_bipartite(G::SimpleGraph)
    
    ADD LATER
    """
    function is_bipartite(G::SimpleGraph)
        color = Dict{Int64, Bool}()
        # for # ADD LATER
        # end
    end
    
    
    #- FUNCTION: `is_cograph`
    """
        is_cograph(G::SimpleGraph)
    
    ADD LATER
    """
    function is_cograph(G::SimpleGraph)
        return !nodesubgraph_is_isomorphic(G, path_graph(4))
    end
end