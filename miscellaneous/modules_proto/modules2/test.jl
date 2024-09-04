using BenchmarkTools
using LinearAlgebra: kron

include("../../modules/KOrthogonalizability.jl")
include("KOrthogonalizability.jl")
is_k_orthogonalizable = KOrthogonalizability.is_k_orthogonalizable
is_k_orthogonalizable2 = KOrthogonalizability2.is_k_orthogonalizable

Y = [-1  -2  -1;
     -1   1  -1;
      1  -1   2];
# Z = [0+0im    2+6im   -4+14im;
#      0+0im   -1+11im   6+6im;
#      4-10im   0+0im    1+15im];
# Z = [-5-2im  0+0im    2+6im   -4+14im;
#      -3-1im  0+0im   -1+11im   6+6im;
#      -2+2im  4-10im   0+0im    1+15im];
Z = [-1.864046107   1.21905858    2.343251229   0.955248594;
      1.118427664  -2.438117161  -7.029753687   0.955248594;
     -1.864046107  -8.533410062   4.686502458  -0.955248594;
     -2.236855328   4.876234321  -9.373004916  -1.432872891];

T = kron(Y, Z);

"""
```julia
julia> is_k_orthogonalizable(T, 5)
(true, [2, 1, 3, 8, 9, 7, 5, 6, 4])
```

```julia
julia> is_k_orthogonalizable2(T, 5)
(true, [1, 2, 3, 7, 9, 8, 4, 6, 5])
"""