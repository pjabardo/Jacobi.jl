""" 
    legendre(x, b)

Compute Legendre polynomials.

Legendre polynomials are Jacobi polynomials with parameters a and b equals to zero.

Same interface as `jacobi`. `dlegendre` computes the derivative of the Legendre 
polynomial.

# Example
```julia-repl
julia> x = legendre(0.4, 5)
0.27063999999999994

julia> x = legendre(2//5, 5)
3383//12500

julia> x = dlegendre(0.4, 5)
-1.3170000000000002

julia> x = dlegendre(2//5, 5)
-1317//1000
```
"""
function legendre(x, n)


    if n==0
        return one(x)
    elseif n==1
        return(x)
    end

    p0 = one(x)
    p1 = x

    for i = 2:n
        p2 = ( (2i-1)*x*p1 - (i-1)*p0 ) / i
        p0 = p1
        p1 = p2
    end

    return p1
end



dlegendre(x, n) = djacobi(x, n)


legendre_zeros(m) = jacobi_zeros(m)

@doc (@doc legendre) dlegendre
@doc (@doc jacobi_zeros!) legendre_zeros



