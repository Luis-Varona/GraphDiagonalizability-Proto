#- Import custom modules for {-1,0,1}-diagonalizability and bandwidth
include("../modules/ZeroOneNegBandwidths.jl")
using .ZeroOneNegBandwidths: is_DiagGraph
DiagGraph = ZeroOneNegBandwidths.GraphObjects.DiagGraph


#- Define broadcasted version of 'is_DiagGraph' to allow for double broadcasting
are_DiagGraphs(input; kwargs...) = is_DiagGraph.(input; kwargs...)
filter_to_diag(input) = map(x -> x[2], filter(x -> x[1], input))


#- Import other Julia modules for loading and saving data
using JLD2: jldsave, load
using PythonCall: Py, pyimport, pylist, pytuple, pyconvert
np = pyimport("numpy")


#-
function DiagGraphs_to_PyArray(graphs::Vector{DiagGraph})::Py
    function DiagGraph_to_PyTuple(g::DiagGraph)::Py
        return pytuple([
            g.graph6_string,
            g.band_zerooneneg,
            g.band_oneneg,
            np.array(g.eigvals),
            np.array(g.eigvecs_zerooneneg),
            np.array(g.eigvecs_oneneg)
        ])
    end
    
    max_order = maximum([g.order for g in graphs])
    string_length = (max_order == 12) ? 12 : (max_order == 13) ? 14 : 17
    PyGraph = np.dtype(pylist([
        ("graph6_string", "<U$string_length"),
        ("band_zerooneneg", "<i8"),
        ("band_oneneg", "O"),
        ("eigvals", "O"),
        ("eigvecs_zerooneneg", "O"),
        ("eigvecs_oneneg", "O")
    ]))
    
    return np.array(pylist(DiagGraph_to_PyTuple.(graphs)), dtype = PyGraph)
end


#- Define main function
function main()
    #- Load Laplacian integral graphs in graph6 format and convert to Laplacians
    graph6_strings_orders12to14 = pyconvert(
        Vector{Vector{String}},
        np.load(
            "graph_data/LaplacianIntegralGraphs/bipartite/" *
            "lap_int_con_bi_graphs_orders12to14.npy",
            allow_pickle = true
        )
    )
    println(length.(graph6_strings_orders12to14))
    
    
    #- Compute {-1,0,1}- and {-1,1}-bandwidths of Laplacian integral graphs
    lap_int_con_bi_graphs_order12 = Vector{Tuple{Bool, DiagGraph}}(
        undef, length(graph6_strings_orders12to14[1])
    )
    for (i, g) in enumerate(graph6_strings_orders12to14[1])
        println("$i (12)")
        lap_int_con_bi_graphs_order12[i] = is_DiagGraph(g)
    end
    diag_con_bi_graphs_order12 = filter_to_diag(lap_int_con_bi_graphs_order12)
    
    
    #-
    lap_int_con_bi_graphs_order13 = Vector{Tuple{Bool, DiagGraph}}(
        undef, length(graph6_strings_orders12to14[2])
    )
    for (i, g) in enumerate(graph6_strings_orders12to14[2])
        println("$i (13)")
        lap_int_con_bi_graphs_order13[i] = is_DiagGraph(
        g, min_zerooneneg = 2, min_oneneg = Inf
    )
    end
    diag_con_bi_graphs_order13 = filter_to_diag(lap_int_con_bi_graphs_order13)
    
    
    #-
    lap_int_con_bi_graphs_order14 = Vector{Tuple{Bool, DiagGraph}}(
        undef, length(graph6_strings_orders12to14[3])
    )
    for (i, g) in enumerate(graph6_strings_orders12to14[3])
        println("$i (14)")
        lap_int_con_bi_graphs_order14[i] = is_DiagGraph(
        g, min_zerooneneg = 2, min_oneneg = 13
    )
    end
    diag_con_bi_graphs_order14 = filter_to_diag(lap_int_con_bi_graphs_order14)
    
    
    #-
    diag_con_bi_graphs_orders12to14 = [
        diag_con_bi_graphs_order12, diag_con_bi_graphs_order13, diag_con_bi_graphs_order14
    ]
    
    
    #-
    if !isdir("graph_data/DiagonalizableGraphs/bipartite")
        mkpath("graph_data/DiagonalizableGraphs/bipartite")
    end
    jldsave(
        "graph_data/DiagonalizableGraphs/bipartite/diag_con_bi_graphs_orders12to14.jld2";
        diag_con_bi_graphs_orders12to14
    )
    
    
    #- Save to named numpy ndarray for Python interoperability
    PyPort_diag_con_bi_graphs_order12 = DiagGraphs_to_PyArray(diag_con_bi_graphs_order12)
    if length(diag_con_bi_graphs_order13) > 0
        PyPort_diag_con_bi_graphs_order13 = DiagGraphs_to_PyArray(
            diag_con_bi_graphs_order13
        )
    else
        PyPort_diag_con_bi_graphs_order13 = np.array(
            pylist(),
            dtype = np.dtype(pylist([
                ("graph6_string", "<U14"),
                ("band_zerooneneg", "<i8"),
                ("band_oneneg", "O"),
                ("eigvals", "O"),
                ("eigvecs_zerooneneg", "O"),
                ("eigvecs_oneneg", "O")
        ])))
    end
    PyPort_diag_con_bi_graphs_order14 = DiagGraphs_to_PyArray(diag_con_bi_graphs_order14)
    
    PyPort_diag_con_bi_graphs_orders12to14 = np.array(
        pylist(
            [PyPort_diag_con_bi_graphs_order12,
            PyPort_diag_con_bi_graphs_order13,
            PyPort_diag_con_bi_graphs_order14]
        ),
        "O"
    )
    np.save(
        "graph_data/DiagonalizableGraphs/bipartite/diag_con_bi_graphs_orders12to14.npy",
        PyPort_diag_con_bi_graphs_orders12to14
    )
end

#- Run main program
# main()