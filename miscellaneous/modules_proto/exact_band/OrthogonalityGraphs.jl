module OrthogonalityGraphs
    #- EXPORTS AND IMPORTS
    export orthograph_complement
    using LinearAlgebra: diagind
    
    
    #- TYPE ALIASES: `ComplexRealInt`, `ComplexRealFloat`
    ComplexRealInt = Union{Complex{Int64}, Int64} # For integer-precision orthogonality
    ComplexRealFloat = Union{Float64, ComplexF64} # For floating-point orthogonality
    
    
    #- FUNCTION: `orthograph_complement`
    """
        orthograph_complement(X::AbstractMatrix{<:ComplexRealInt})::BitMatrix
        orthograph_complement(X::AbstractMatrix{<:ComplexRealFloat})::BitMatrix
    
    Compute the adjacency matrix of the orthogonality graph complement of a matrix.
    
    # Arguments
    - `X::AbstractMatrix{<:ComplexRealInt}`: a (potentially complex) integer matrix. (Exact)
        orthogonality between columns is required.)
    - `X::AbstractMatrix{<:ComplexRealFloat}`: a (potentially complex) floating-point
        matrix. (Only approximate orthogonality is required.)
    
    # Returns
    - `A::BitMatrix`: the adjacency matrix of the orthogonality graph complement of `X`.
    
    # Examples
    Process a `4`-orthogonalizable integer matrix:
    ```jldoctest
    julia> X = [-2  -2   3  -4   1;
                -2  -3  -2   2  -2;
                -6   1   0  -1   1];
    julia> orthograph_complement(X)
    5×5 BitMatrix:
     0  1  1  1  1
     1  0  0  1  1
     1  0  0  1  1
     1  1  1  0  1
     1  1  1  1  0
    ```
    
    Process a quasi-orthogonal (complex) floating-point matrix:
    ```jldoctest
    julia> X = [ 7.1093617+8.2942553im  -3.2708609-1.6354304im   3.8275079-3.8275079im;
                 1.7773404+0.5924468im   0.5451435-9.2674392im        -0.0-11.4825237im;
                -5.9244681+5.9244681im  -0.5451435-1.6354304im  -3.8275079-3.8275079im;
                 -2.962234+2.962234im    1.6354304-8.7222957im  -7.6550158-3.8275079im];
    julia> orthograph_complement(X)
    3×3 BitMatrix:
     0  1  0
     1  0  1
     0  1  0
    ```
    """
    function orthograph_complement(X::AbstractMatrix{<:ComplexRealInt})::BitMatrix
        A = (X' * X .!= 0) # Take the Gram matrix of the columns of `X`
        A[diagind(A)] .= 0 # Remove self-loops from the resulting graph
        return A
    end
    
    function orthograph_complement(X::AbstractMatrix{<:ComplexRealFloat})::BitMatrix
        # Take the Gram matrix of the columns of `X` (with a tolerance for orthogonality)
        A = .!isapprox.(X' * X, 0, atol = 1e-8)
        A[diagind(A)] .= 0 # Remove self-loops from the resulting graph
        return A
    end
end