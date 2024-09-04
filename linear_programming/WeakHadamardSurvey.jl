#-
include("NormalizedWeakHadamards.jl")
include("../modules/ZeroOneNegBandwidths.jl")

using .NormalizedWeakHadamards: normalized_weak_hadamards
using .ZeroOneNegBandwidths: is_DiagGraph
using .ZeroOneNegBandwidths.GraphObjects: DiagGraph, graph6_to_laplacian

#-
using JLD2: save, load
using PythonCall
np = pyimport("numpy")
sc = pyimport("scipy")

#-
function weak_hads()
    norm_weak_hads_orders1to6 = Vector{Vector{Matrix{Int64}}}(undef, 6)
    for n in 1:6
        norm_weak_hads_orders1to6[n] = normalized_weak_hadamards(n)
    end
    
    isdir("linear_programming/data") || mkpath("linear_programming/data")
    save(
        "linear_programming/data/norm_weak_hads_orders1to6.jld2",
        "norm_weak_hads_orders1to6", norm_weak_hads_orders1to6)
    
    PyPort_norm_weak_hads_orders1to6 = np.array(
        [np.array(pylist([np.array(i) for i in j]))
        for j in norm_weak_hads_orders1to6])
    MatPort_norm_weak_hads_orders1to6 = np.array([
        bundle.swapaxes(1, 2).swapaxes(0, 2)
        for bundle in PyPort_norm_weak_hads_orders1to6])
    sc.io.savemat(
        "linear_programming/data/norm_weak_hads_orders1to6.mat",
        Dict(
            [(
                "norm_weak_hads_order$n",
                MatPort_norm_weak_hads_orders1to6[n - 1])
             for n in 1:6]))
end

#-
function laplacian_to_adjacency(L::AbstractMatrix{Int64})
    return Diagonal(diag(L)) - L
end

#-
function adjacencies()
    graph6_strings_orders3to6 = pyconvert(
        Vector{Vector{String}},
        np.load(
            "linear_programming/data/graph6_strings_orders3to6.npy",
            allow_pickle = true))
    laplacians_orders3to6 = map.(
        graph6_to_laplacian, graph6_strings_orders3to6)
    
    non_WHD_adjacencies_orders3to6 = map.(
        laplacian_to_adjacency,
        filter.(
            x -> oneneg_bandwidth(x, include_zero = true)[1] > 2,
            laplacians_orders3to6))
    save(
        "linear_programming/data/non_WHD_adjacencies_orders3to6.jld2",
        "non_WHD_adjacencies_orders3to6", non_WHD_adjacencies_orders3to6)
    
    PyPort_non_WHD_adjacencies_orders3to6 = np.array(
        [np.array(pylist([np.array(i) for i in j]))
        for j in non_WHD_adjacencies_orders3to6])
    MatPort_non_WHD_adjacencies_orders3to6 = np.array([
        bundle.swapaxes(1, 2).swapaxes(0, 2)
        for bundle in PyPort_non_WHD_adjacencies_orders3to6])
    sc.io.savemat(
        "linear_programming/data/non_WHD_adjacencies_orders3to6.mat",
        Dict(
            [(
                "non_WHD_adjacencies_order$n",
                MatPort_non_WHD_adjacencies_orders3to6[n - 3])
             for n in 3:6]))
end

#-
function main()
    weak_hads()
    adjacencies()
end

main()