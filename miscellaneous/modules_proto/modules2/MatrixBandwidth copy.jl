module MatrixBandwidth
    export saxe_search
    
    using Combinatorics: permutations
    using DataStructures: Queue
    using LinearAlgebra: diagind
    
    Base.popfirst!(s::Queue) = popfirst!(s.store)
    function Base.push!(q::Queue, x)
        push!(q.store, x)
        return q
    end
    
    function saxe_search(A::BitMatrix, k::Int64)
        if any(sum(A, dims = 2) .> 2k - 2)
            return (false, Int64[])
        end
        
        n = size(A, 2)
        
        R = collect(Iterators.flatten([permutations(1:n, i) for i in 0:n]))
        Q = Queue{Int64}()
        arr = [[false, Int64[]] for _ in 1:length(R)]
        
        push!(Q, 1)
        arr[1] = [true, collect(1:n)]
        
        while !isempty(Q)
            r = popfirst!(Q)
            unplaced_nodes = arr[r][2]
            # if R[r][1:min(2, end)] == [2,1]
            #     println(R[r])
            #     println(unplaced_nodes)
            # end
            
            for node in unplaced_nodes
                region = vcat(R[r], node)
                s = findfirst(x -> x == region, R)
                
                if !isnothing(s) && !arr[s][1]
                    arr[s][1] = true
                    arr[s][2] = filter(x -> x != node, arr[r][2])
                    
                    if !any(A[node, region[1:max(end - k, 0)]])
                        if isempty(arr[s][2])
                            order = R[s]
                            if order[1] > order[end]
                                order = reverse(order) # To better preserve the original ordering
                            end
                            return (true, order)
                        end
                        
                        push!(Q, s)
                    end
                end
            end
        end
        
        return (false, Int64[])
    end
end