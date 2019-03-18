using Polynomials


"""
Calculate Chebyshev polynomials coefficients

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
function poly_chebyshev(n, ::Type{T}=Int, var=:x) where {T<:Number}

    if n==0
        return Poly{T}([one(T)], var)
    elseif n==1
        return Poly{T}([zero(T), one(T)], var) 
    end

    px = Poly{T}([zero(T), one(T)], var)
    T0 = Poly{T}([one(T)], var)
    T1 = px
    
    for i = 2:n
        T2 = 2*px*T1 - T0
        T0 = T1
        T1 = T2
    end

    return T1

end

"""
Calculate derivative of Chebyshev polynomials coefficients

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
poly_dchebyshev(n, ::Type{T}=Int, var=:x) where {T<:Number} = polyder(poly_chebyshev(n,T, var))



"""
Calculate the coefficients of Chebyshev polynomials of the second kind

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
function poly_chebyshev2(n, ::Type{T}=Int, var=:x) where {T<:Number}

    if n==0
        return Poly{T}([one(T)], var)
    elseif n==1
        return Poly{T}([zero(T), 2*one(T)], var) 
    end

    px = Poly{T}([zero(T), one(T)], var)
    U0 = Poly{T}([one(T)], var)
    U1 = Poly{T}([zero(T), 2*one(T)], var)
    
    for i = 2:n
        U2 = 2*px*U1 - U0
        U0 = U1
        U1 = U2
    end

    return U1

end

"""
Calculate the coefficients of the derivative of Chebyshev polynomials of the second kind

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
poly_dchebyshev2(n, ::Type{T}=Int, var=:x) where {T<:Number} = polyder(poly_chebyshev2(n,T, var))


"""
Calculate the coefficients of Legendre polynomials 

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
function poly_legendre(n, ::Type{T}=Float64, var=:x) where {T<:Number}


    if n==0
        return Poly{T}([one(T)], var)
    elseif n==1
        return Poly{T}([zero(T), one(T)], var) 
    end

    px = Poly{T}([zero(T), one(T)], var)
    p0 = Poly{T}([one(T)], var)
    p1 = px

    for i = 2:n
        p2 = ( (2i-1)*px*p1 - (i-1)*p0 ) / i
        p0 = p1
        p1 = p2
    end

    return p1
end


"""
Calculate the coefficients of the derivativ of Legendre polynomials 

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
poly_dlegendre(n, ::Type{T}=Int, var=:x) where {T<:Number} = polyder(poly_legendre(n,T, var))





"""
Calculate the coefficients of Jacobi polynomials 

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
function poly_jacobi(n, a, b, ::Type{T}=Float64, var=:x) where {T<:Number}

    ox = one(T)
    zx = zero(T)
    if n==0
        return Poly{T}([one(T)], var)
    elseif n==1
        return Poly{T}([(a-b)/2, (a+b+2)/2], var)
    end
    p0 = Poly{T}([one(T)], var)
    p1 = Poly{T}([(a-b)/2, (a+b+2)/2], var)
    px = Poly{T}([zero(T), one(T)], var)
    
    for i = 1:(n-1)
	a1 = 2*(i+1)*(i+a+b+1)*(2*i+a+b);
	a2 = (2*i+a+b+1)*(a*a-b*b);
	a3 = (2*i+a+b)*(2*i+a+b+1)*(2*i+a+b+2);
	a4 = 2*(i+a)*(i+b)*(2*i+a+b+2);
	p2 = ox/a1*( (a2 + a3*px)*p1 - a4*p0);

        p0 = p1
        p1 = p2
    end

    return p1
end


"""
Calculate the coefficients of the derivative of Jacobi polynomials 

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial
"""
poly_djacobi(n, a, b, ::Type{T}=Float64, var=:x) where {T<:Number} = polyder(poly_jacobi(n,a,b,T,var))
