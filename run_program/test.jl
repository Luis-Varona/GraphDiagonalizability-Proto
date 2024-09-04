"""
Computes the {-1,0,1}- and {-1,1}-bandwidths of the Laplacian matrix in
"example_matrix.txt". To test your own matrix, uncomment line 16, comment out
line 17, paste your desired matrix into "your_matrix.txt", and run the program.
(However, make sure your working directory is the entire project directory
and not the "run_program" subdirectory, or you may run into problems.)
"""
#- Import required modules
using DelimitedFiles: readdlm # To read the Laplacian matrix from a text file
include("../modules/ZeroOneNegBandwidths.jl")
using .ZeroOneNegBandwidths: is_DiagGraph # To compute S-bandwidth

#-
function main()
    # Load the Laplacian matrix from a text file
    # L = readdlm("run_program/your_matrix.txt", Int64) # To actually test (uncomment this)
    L = readdlm("run_program/example_matrix.txt", Int64) # The example (comment this out)
    
    # Replace these if you already have some lower bound on the S-bandwidths
    min_zerooneneg = 1
    min_oneneg = 1
    
    # Identifies the {-1,0,1}- and {-1,1}-spectra of the Laplacian matrix
    (bool, data) = is_DiagGraph(L, min_zerooneneg = min_zerooneneg, min_oneneg = min_oneneg)
    
    # Display the Laplacian matrix
    println("Laplacian matrix:")
    display(L)
    
    # Display the {-1,0,1}- and {-1,1}-bandwidths
    println("\n{-1,0,1}-bandwidth: $(data.band_zerooneneg)")
    println("{-1,1}-bandwidth: $(data.band_oneneg)\n")
    
    # Display the Laplacian spectra if it is a {-1,0,1}-diagonalizable graph
    if bool
        # Display the eigenvalues
        println("Eigenvalues:")
        display(data.eigvals)
        
        # Display the associated {-1,0,1}-eigenvectors
        println("\nMatrix of {-1,0,1}-eigenvectors:")
        display(data.eigvecs_zerooneneg)
        
        # Display the associated {-1,1}-eigenvectors if applicable
        if data.band_oneneg < Inf
            println("\nMatrix of {-1,1}-eigenvectors:")
            display(data.eigvecs_oneneg)
        end
    end
end

main()