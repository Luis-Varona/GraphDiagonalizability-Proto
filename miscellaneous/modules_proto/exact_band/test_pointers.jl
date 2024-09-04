function test_pointers()
    a = [1, 2, 3, 4, 5]
    b = a
    return (a, b)
end

function main()
    (a, b) = test_pointers()
    a[2] = 10
    println(a)
    println(b)
end

main()