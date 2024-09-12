#- Imports
include("UniqueWithTolerance.jl")
using .UniqueWithTolerance: uniquetol
using Distributions: Uniform
using Random: shuffle!

function main()
    #-
    numOnes = rand(1:13)
    ctsInitial = Vector{Int64}(undef, 100)
    ctsInitial[1:numOnes] .= 1
    ctsInitial[(numOnes + 1):100] .= rand(2:13, 100 - numOnes)

    #-
    rangeVals = (rand(-87:-13), rand(13:87))
    valsInitial = rand(Uniform(rangeVals[1], rangeVals[2]), 100)

    #-
    rangeShifts = (rand(-334:-87), rand(87:334))
    shifts = rand(rangeShifts[1]:rangeShifts[2], 100)
    valsInitial .+= shifts

    #-
    diffs = map(x -> rand(Uniform(-1e-7, 1e-7), x), ctsInitial)
    testArr = vcat([v .+ ct for (v, ct) in zip(valsInitial, diffs)]...)
    shuffle!(testArr)

    #-
    (u, ix, cts) = uniquetol(testArr; return_indices=true, return_counts=true)
    spaces = diff(u)
    
    #-
    return (testArr=testArr, u=u, ix=ix, cts=cts, spaces=spaces)
end

(testArr, u, ix, cts, spaces) = main()