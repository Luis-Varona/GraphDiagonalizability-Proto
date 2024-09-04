function result = ullmanSearch(motif, host, assignments)
    l = length(assignments);
    edges = motif.Edges;
    
    for k = 1:size(edges, 1)
        edge = edges(k, 1);
        node1 = edge.EndNodes(1);
        node2 = edge.EndNodes(2);
        if node1 < l && node2 < l
            hostEdge = sort(assignments([node1, node2]));
            if not(ismember(hostEdge, host.Edges.EndNodes, 'rows'))
                result = false;
                return;
            end
        end
    end
    
    if k == size(edges, 1)
        result = true;
        return;
    end
    
    meow
end