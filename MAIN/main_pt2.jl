#- Import modules to test for {-1,0,1}-diagonalizability and bandwidth
include("../modules/ZeroOneNegBandwidths.jl")
using .ZeroOneNegBandwidths: is_DiagGraph
using .ZeroOneNegBandwidths.GraphObjects: DiagGraph, DiagGraphs_to_PyArray


#- Define additional functions
# A broadcasted version of 'is_DiagGraph' (to allow for double broadcasting)
are_DiagGraphs(input; kwargs...) = is_DiagGraph.(input; kwargs...)

# Filters an array of 'is_DiagGraph' results to {-1,0,1}-diagonalizable graphs only
filter_to_diag(input) = map(x -> x[2], filter(x -> x[1], input))


#- Import other Julia modules for loading and saving data
using JLD2
using PythonCall
np = pyimport("numpy")


#- Define main function
function main()
    #- Load all Laplacian integral simple connected graphs up to order 11 in graph6 format
    graph6_strings_orders1to10 = pyconvert(
        Vector{Vector{String}},
        np.load(
            "graph_data/LaplacianIntegralGraphs/lap_int_con_graphs_orders1to10.npy",
            allow_pickle = true
        )
    )
    graph6_strings_order11 = pyconvert(
        Vector{String},
        np.load("graph_data/LaplacianIntegralGraphs/lap_int_con_graphs_order11.npy")
    )
    graph6_strings_orders1to11 = vcat(graph6_strings_orders1to10, [graph6_strings_order11])
    
    
    #- Compute the S-bandwidths of all Laplacian integral graphs up to order 8
    lap_int_con_graphs_orders1to8 = are_DiagGraphs.(graph6_strings_orders1to11[1:8])
    diag_con_graphs_orders1to8 = filter_to_diag.(lap_int_con_graphs_orders1to8)
    
    
    #= Compute the S-bandwidths of all Laplacian integral graphs of orders 9, 10, and 11.
    Theoretical lower bounds from Johnston & Plosker 2024 are included to save time. =#
    lap_int_con_graphs_order9 = are_DiagGraphs(
        graph6_strings_orders1to11[9], min_zerooneneg = 2, min_oneneg = Inf
    )
    diag_con_graphs_order9 = filter_to_diag(lap_int_con_graphs_order9)
    
    # Iterate over orders 10 and 11 instead of broadcasting to save RAM
    lap_int_con_graphs_order10 = Vector{Tuple{Bool, DiagGraph}}(
        undef, length(graph6_strings_orders1to11[10])
    )
    for (i, g) in enumerate(graph6_strings_orders1to11[10])
        println("$i (10)") ###
        lap_int_con_graphs_order10[i] = is_DiagGraph(g, min_zerooneneg = 2, min_oneneg = 9)
    end
    diag_con_graphs_order10 = filter_to_diag(lap_int_con_graphs_order10)
    
    lap_int_con_graphs_order11 = Vector{Tuple{Bool, DiagGraph}}(
        undef, length(graph6_strings_orders1to11[11])
    )
    for (i, g) in enumerate(graph6_strings_orders1to11[11])
        println("$i (11)") ###
        lap_int_con_graphs_order11[i] = is_DiagGraph(
            g, min_zerooneneg = 2, min_oneneg = Inf
        )
    end
    diag_con_graphs_order11 = filter_to_diag(lap_int_con_graphs_order11)
    
    
    #- Combine all {-1,0,1}-diagonalizable graphs up to order 11 into a single array
    diag_con_graphs_orders1to11 = vcat(
        diag_con_graphs_orders1to8,
        [diag_con_graphs_order9, diag_con_graphs_order10, diag_con_graphs_order11]
    )
    
    
    #- Save results to a JLD2 file
    if !isdir("graph_data/DiagonalizableGraphs")
        mkpath("graph_data/DiagonalizableGraphs")
    end
    jldsave(
        "graph_data/DiagonalizableGraphs/diag_con_graphs_orders1to11.jld2";
        diag_con_graphs_orders1to11
    )
    
    #- Save results to a named numpy ndarray for Python interoperability
    PyPort_diag_con_graphs_orders1to11 = DiagGraphs_to_PyArray.(diag_con_graphs_orders1to11)
    np.save(
        "graph_data/DiagonalizableGraphs/diag_con_graphs_orders1to11.npy",
        PyPort_diag_con_graphs_orders1to11
    )
end


#- Run the main program
main()