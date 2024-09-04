module ProtoEigvals
    #- EXPORTS AND IMPORTS
    export algorithm1A, algorithm1B, algorithm2
    
    using LinearAlgebra: Symmetric, eigvals
    using DataStructures: Accumulator, counter
    using FreqTables: freqtable
    
    
    #- FUNCTION: `algorithm1A`
    function algorithm1A(
        L::AbstractMatrix{Int64}
    )::Tuple{Bool, Vector{Int64}, Vector{Tuple{Int64, Int64}}}
        
        # Computing eigenvalues is more efficient for `Symmetric`-type matrices
        Λ_float = eigvals(Symmetric(L, :L))
        Λ = round.(real.(Λ_float))
        is_integral = isapprox(Λ, Λ_float)
        
        if !is_integral
            Λ = Int64[]
            M = Tuple{Int64, Int64}[]
        else
            # Count all multiplicitiess and sort unique eigenvalues by multiplicity
            λ_table = sort(freqtable(Int64.(Λ)))
            M = [(λ, μ) for (λ, μ) in zip(names(λ_table)[1], collect(λ_table))]
            
            # Sort all eigenvalues, including repetitions, by multiplicity
            Λ = vcat(fill.(first.(M), last.(M))...)
        end
        
        return (is_integral, Λ, M)
    end
    
    
    #- FUNCTION: `algorithm1B`
    function algorithm1B(
        L::AbstractMatrix{Int64}
    )::Tuple{Bool, Vector{Int64}, Vector{Pair{Int64, Int64}}}
        
        # Computing eigenvalues is more efficient for `Symmetric`-type matrices
        Λ_float = eigvals(Symmetric(L, :L))
        Λ = round.(real.(Λ_float))
        is_integral = isapprox(Λ, Λ_float)
        
        if !is_integral
            Λ = Int64[]
            M = Pair{Int64, Int64}[]
        else
            # Count all multiplicitiess and sort unique eigenvalues by multiplicity
            λ_table = sort(freqtable(Int64.(Λ)))
            M = [λ => μ for (λ, μ) in zip(names(λ_table)[1], collect(λ_table))]
            
            # Sort all eigenvalues, including repetitions, by multiplicity
            Λ = vcat(fill.(first.(M), last.(M))...)
        end
        
        return (is_integral, Λ, M)
    end
    
    
    #- FUNCTION: `algorithm2`
    function algorithm2(
        L::AbstractMatrix{Int64}
    )::Tuple{Bool, Vector{Int64}, Vector{Pair{Int64, Int64}}}
        
        # Computing eigenvalues is more efficient for `Symmetric`-type matrices
        Λ_float = eigvals(Symmetric(L, :L))
        Λ = round.(real.(Λ_float))
        is_integral = isapprox(Λ, Λ_float)
        
        if !is_integral
            Λ = Int64[]
            λ_counter = Pair{Int64, Int64}[]
        else
            # Count all multiplicitiess and sort unique eigenvalues by multiplicity
            λ_counter = sort_counter(counter(Int64.(Λ)))
                        
            # Sort all eigenvalues, including repetitions, by multiplicity
            Λ = vcat(fill.(first.(λ_counter), last.(λ_counter))...)
        end
        
        return (is_integral, Λ, λ_counter)
    end
    
    
    #- Misc
    # most_common(c::Accumulator) = most_common(c, length(c))
    # most_common(c::Accumulator, k) = sort(collect(c), by = kv -> kv[2], rev = true)[1:k]
    sort_counter(c::Accumulator; rev::Bool = false) = sort(collect(c), by = kv -> kv[2], rev = rev)
end

# using DataStructures ###
# using FreqTables ###
# const arr = rand(-5:5, 17) ###

using BenchmarkTools
using .ProtoEigvals

const L = [[15 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 15 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 14  0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1  0 14 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 11  0  0  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0 11  0  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0 11  0  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0  0 11  0 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1  0  0  0  0 11 -1 -1 -1 -1 -1 -1 -1]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  9  0  0  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  9  0  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  9  0  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  9  0  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  9  0  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  9  0]
    [-1 -1 -1 -1 -1 -1 -1 -1 -1  0  0  0  0  0  0  9]];