""" 
Compute Legendre polynomials

Legendre polynomials are Jacobi polynomials with parameters a and b equals to zero.

Same interface as `jacobi`. `dlegendre` computes the derivative of the Legendre 
polynomial.
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



