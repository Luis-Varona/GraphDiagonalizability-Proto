#-
include("GraphObjects.jl")
using .GraphObjects: DiagGraph
using Intervals: (..)


#-
function main()
    println("NOTE: This is not a real graph, just a test case of the `DiagGraph` type.")
    Γ = DiagGraph(
        [3 -2 -1; -2 2 0; -1 0 1],
        1 .. 3,
        Inf,
        [0, 1, 2],
        [1 1 1; 1 0 -1; 1 -1 0],
        [missing;;],
    )
    
    
    println("\n\nREPR")
    display(Γ)
    
    
    println("\n\nFIELDS")
    print("Laplacian matrix: "); display(Γ.laplacian_matrix)
    print("\n{-1,0,1}-bandwidth: ");display(Γ.band_zerooneneg)
    print("{-1,1}-bandwidth: "); display(Γ.band_oneneg)
    print("\nEigenvalues: "); display(Γ.eigvals)
    print("\n{-1,0,1}-eigenvectors: "); display(Γ.eigvecs_zerooneneg)
    print("\n{-1,1}-eigenvectors: "); display(Γ.eigvecs_oneneg)
    
    
    println("\n\nMETHODS (Part 1)")
    print("Adjacency matrix: "); display(Γ.adjacency_matrix)
    print("\ngraph6 string: "); println("[Feature not yet implemented]")
    
    println("\n\nMETHODS (Part 2)")
    print("Order: "); display(Γ.order)
    print("Size: "); display(Γ.size)
    print("Density: "); display(Γ.density)
    print("Weighted density: "); display(Γ.weighted_density)
    print("Average degree: "); display(Γ.average_degree)
    print("Average weighted degree: "); display(Γ.average_weighted_degree)
    
    println("⤱\nMETHODS (Part 3)")
    print("Is weighted: "); display(Γ.is_weighted)
    print("Is negatively weighted: "); display(Γ.is_negatively_weighted)
    print("Is connected: "); display(Γ.is_connected)
    print("Is regular: "); display(Γ.is_regular)
    print("Is weighted regular: "); display(Γ.is_weighted_regular)
    print("Is bipartite: "); display(Γ.is_bipartite)
    
    println("\nMETHODS (Part 4)")
    print("Is cograph: "); display(Γ.is_cograph)
    print("Is Cartesian product: "); println("[Feature not yet implemented]")
    print("Is Cartesian product complement: ") ; println("[Feature not yet implemented]")
end


#-
main()