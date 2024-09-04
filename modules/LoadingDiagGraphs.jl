module LoadingDiagGraphs
    #- EXPORTS AND IMPORTS
    export load_file
    
    using JLD2: jldopen
    
    include("GraphObjects.jl")
    using .GraphObjects: DiagGraph
    
    #- FUNCTION: `load_file`
    function load_file(path::String)
        return jldopen(
            path, "r";
            typemap = Dict(
                "Main.ZeroOneNegBandwidths.GraphObjects.DiagGraph" => GraphObjects.DiagGraph
            )
        )
    end
end