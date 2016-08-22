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

function legendre!{T<:Number}(x::AbstractArray{T}, n, y::AbstractArray{T})

    m = length(x)
    for i = 1:m
        y[i] = legendre(x[i], n)
    end
    return y
end

legendre{T<:Number}(x::AbstractArray{T}, n) = legendre!(x, n, zeros(x))


dlegendre(x, n) = djacobi(x, n)


function dlegendre!{T<:Number}(x::AbstractArray{T}, n, y::AbstractArray{T})

    m = length(x)
    for i = 1:m
        y[i] = dlegendre(x[i], n)
    end
    return y
end

dlegendre{T<:Number}(x::AbstractArray{T}, n) = dlegendre!(x, n, zeros(x))
legendre_zeros(m) = jacobi_zeros(m)

@doc (@doc legendre) legendre!
@doc (@doc legendre) dlegendre
@doc (@doc legendre) dlegendre!
@doc (@doc jacobi_zeros!) legendre_zeros



