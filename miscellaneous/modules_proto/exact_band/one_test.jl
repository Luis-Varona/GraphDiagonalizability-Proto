#-
using Graphs: SimpleGraph
using MolecularGraph: nodesubgraph_isomorphisms
include("LaplacianSpectra.jl")

#-
const n = 12
E = Vector{Vector{Float64}}(collect(LaplacianSpectra.laplacian_eigvecs(n)))
F = filter(x -> !any(x .== 0), E)
E = hcat(E...)
F = hcat(F...)

A = .!isapprox.(E' * E, 0, atol = 1e-8)
B = .!isapprox.(F' * F, 0, atol = 1e-8)

H = SimpleGraph(zeros(Bool, n - 1, n - 1))

#-
# G1 = SimpleGraph(A)
# M1 = nodesubgraph_isomorphisms(G1, H)
# m1 = iterate(M1)
# if !isnothing(m1)
#     m1 = first.(sort(collect(m1[1])))
#     Q1 = E[:, m1]
# else
#     Q1 = Matrix{Int64}(undef, n, 0)
# end

#-
G2 = SimpleGraph(B)
M2 = nodesubgraph_isomorphisms(G2, H)
m2 = iterate(M2)
if !isnothing(m2)
    m2 = first.(sort(collect(m2[1])))
    Q2 = F[:, m2]
else
    Q2 = Matrix{Int64}(undef, n, 0)
end
