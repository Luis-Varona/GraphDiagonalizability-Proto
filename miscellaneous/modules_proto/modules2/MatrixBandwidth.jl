module MatrixBandwidth
    export saxe_search
    
    using Combinatorics: permutations
    using DataStructures: Queue
    # using LinearAlgebra: diagind
    
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
        # B = zeros(Bool, n, n)
        # for i in vcat((1 - k):-1, 1:(k - 1)) # Make all possible columns non-orthogonal
        #     B[diagind(B, i)] .= true
        # end
        
        R = collect(Iterators.flatten([permutations(1:n, i) for i in 0:(k - 1)]))
        Q = Queue{Int64}()
        arr = [[false, Int64[]] for _ in 1:length(R)]
        
        push!(Q, 1)
        arr[1] = [true, collect(1:n)]
        
        while !isempty(Q)
            r = popfirst!(Q)
            unplaced_nodes = arr[r][2]
            # println(unplaced_nodes) ###
            
            # f = R[r] ###
            # if !isempty(f)
            #     println(A[unplaced_nodes, f])
            #     filter!(x -> any(A[x, f]), unplaced_nodes) ###
            # end
            # println(unplaced_nodes) ###
            
            for node in unplaced_nodes
                region = vcat(R[r], node)
                s = findfirst(x -> x == region, R)
                
                if !isnothing(s) && !arr[s][1]
                    arr[s][1] = true
                    # slice = 1:length(region)
                    # println(slice)
                    # filter!(x -> x != node, arr[r][2])
                    # arr[s][2] = copy(arr[r][2])
                    
                    # if !any(A[node, region[1:max(end - k, 0)]]) # Maybe the problem?
                    # display(B[slice, slice] - A[region, region])
                    # if all(B[slice, slice] - A[region, region] .>= 0)
                        # arr[s][1] = true
                    
                    filter!(x -> x != node, arr[r][2])
                    arr[s][2] = copy(arr[r][2])
                    # arr[s][2] = filter(x -> x != node, arr[r][2])
                    println(arr[s][2])
                    
                    if isempty(arr[s][2])
                        order = R[s]
                        if order[1] > order[end]
                            order = reverse(order) # To better preserve the original ordering
                        end
                        return (true, order)
                    end
                    
                    push!(Q, s)
                    # end
                end
            end
        end
        
        return (false, Int64[])
    end
end