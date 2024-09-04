#-
motzkin_humps_og(n::Int64)::Int64 = sum([binomial(n, k)*binomial(n - k, k)/2 for k in 1:n])
motzkin_humps(n::Int64)::Int64 = sum([binomial(n, k)*binomial(n - k, k)/2 for k in 1:floor(Int64, n/2)])

#-
function test(n::Int64)::Int64
    inner_loop(r::Int64) = sum(
        [factorial(r) / (factorial(i) * factorial(r - 2i + 1) * factorial(i - 1)) for i in 1:ceil(Int64, r/2)]
    )
    return Int64(sum([inner_loop(k) for k in 1:(n - 1)]))
end