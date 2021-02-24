using Polynomials


"""
    poly_chebyshev(n, [::Type{T}=Int, [var=:x]])

Create a Chebyshev polynomial of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_chebyshev(4)
Polynomials.Polynomial(1 - 8*x^2 + 8*x^4)

julia> p2 = poly_chebyshev(4, BigInt)
Polynomials.Polynomial(1 - 8*x^2 + 8*x^4)

julia> p3 = poly_chebyshev(4, Float64, :y)
Polynomials.Polynomial(1.0 - 8.0*y^2 + 8.0*y^4)

```

"""
function poly_chebyshev(n, ::Type{T}=Int, var=:x) where {T<:Number}

    if n==0
        return Polynomial{T}([one(T)], var)
    elseif n==1
        return Polynomial{T}([zero(T), one(T)], var) 
    end

    px = Polynomial{T}([zero(T), one(T)], var)
    T0 = Polynomial{T}([one(T)], var)
    T1 = px
    
    for i = 2:n
        T2 = 2*px*T1 - T0
        T0 = T1
        T1 = T2
    end

    return T1

end

"""
    poly_dchebyshev(n, [::Type{T}=Int, [var=:x]])

Create the derivative of Chebyshev polynomial of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_dchebyshev(4)
Polynomials.Polynomial(-16*x + 32*x^3)

julia> p2 = poly_dchebyshev(4, BigInt)
Polynomials.Polynomial(-16*x + 32*x^3)

julia> p3 = poly_dchebyshev(4, Float64, :y)
Polynomials.Polynomial(-16.0*y + 32.0*y^3)
```
"""
poly_dchebyshev(n, ::Type{T}=Int, var=:x) where {T<:Number} = derivative(poly_chebyshev(n,T, var))



"""
    poly_chebyshev2(n, [::Type{T}=Int, [var=:x]])

Create a Chebyshev of the second kind  of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_chebyshev2(4)
Polynomials.Polynomial(1 - 12*x^2 + 16*x^4)

julia> p2 = poly_chebyshev2(4, BigInt)
Polynomials.Polynomial(1 - 12*x^2 + 16*x^4)

julia> p3 = poly_chebyshev2(4, Float64, :y)
Polynomials.Polynomial(1.0 - 12.0*y^2 + 16.0*y^4)
```

"""
function poly_chebyshev2(n, ::Type{T}=Int, var=:x) where {T<:Number}

    if n==0
        return Polynomial{T}([one(T)], var)
    elseif n==1
        return Polynomial{T}([zero(T), 2*one(T)], var) 
    end

    px = Polynomial{T}([zero(T), one(T)], var)
    U0 = Polynomial{T}([one(T)], var)
    U1 = Polynomial{T}([zero(T), 2*one(T)], var)
    
    for i = 2:n
        U2 = 2*px*U1 - U0
        U0 = U1
        U1 = U2
    end

    return U1

end

"""
    poly_dchebyshev2(n, [::Type{T}=Int, [var=:x]])

Create a derivative of Chebyshev polynomial of the second kind  of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_dchebyshev2(4)
Polynomials.Polynomial(-24*x + 64*x^3)

julia> p2 = poly_dchebyshev2(4, BigInt)
Polynomials.Polynomial(-24*x + 64*x^3)

julia> p3 = poly_dchebyshev2(4, Float64, :y)
Polynomials.Polynomial(-24.0*y + 64.0*y^3)
```

"""
poly_dchebyshev2(n, ::Type{T}=Int, var=:x) where {T<:Number} = derivative(poly_chebyshev2(n,T, var))


"""
    poly_legendre(n, [::Type{T}=Int, [var=:x]])

Create a Legendre polynomial of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_legendre(4)
Polynomials.Polynomial(0.375 - 3.75*x^2 + 4.375*x^4)

julia> p2 = poly_legendre(4, Rational{Int})
Polynomials.Polynomial(3//8 - 15//4*x^2 + 35//8*x^4)

julia> p3 = poly_legendre(4, BigFloat, :y)
Polynomials.Polynomial(0.375 - 3.75*y^2 + 4.375*y^4)
```

"""
function poly_legendre(n, ::Type{T}=Float64, var=:x) where {T<:Number}


    if n==0
        return Polynomial{T}([one(T)], var)
    elseif n==1
        return Polynomial{T}([zero(T), one(T)], var) 
    end

    px = Polynomial{T}([zero(T), one(T)], var)
    p0 = Polynomial{T}([one(T)], var)
    p1 = px

    for i = 2:n
        p2 = ( (2i-1)*px*p1 - (i-1)*p0 ) / i
        p0 = p1
        p1 = p2
    end

    return p1
end


"""
    poly_dlegendre(n, [::Type{T}=Int, [var=:x]])

Create a derivative of Legendre polynomial of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_dlegendre(4)
Polynomials.Polynomial(-7.5*x + 17.5*x^3)

julia> p2 = poly_dlegendre(4, Rational{Int})
Polynomials.Polynomial(-15//2*x + 35//2*x^3)

julia> p3 = poly_dlegendre(4, BigFloat, :y)
Polynomials.Polynomial(-7.5*y + 17.5*y^3)
```

"""
poly_dlegendre(n, ::Type{T}=Int, var=:x) where {T<:Number} = derivative(poly_legendre(n,T, var))





"""
    poly_jacobi(n, [::Type{T}=Int, [var=:x]])

Create a Jacobi polynomial of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `a` Weight α of the Jacobi polynomial
 * `b` Weight β of the Jacobi polynomial
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_jacobi(4, 1, 1)
Polynomials.Polynomial(0.625 - 8.75*x^2 + 13.125*x^4)

julia> p2 = poly_jacobi(4, 1, 1, Rational{Int})
Polynomials.Polynomial(5//8 - 35//4*x^2 + 105//8*x^4)

julia> p3 = poly_jacobi(4, 1, 1, BigFloat, :y)
Polynomials.Polynomial(0.625 - 8.75*y^2 + 13.125*y^4)
```
"""
function poly_jacobi(n, a, b, ::Type{T}=Float64, var=:x) where {T<:Number}

    ox = one(T)
    zx = zero(T)
    if n==0
        return Polynomial{T}([one(T)], var)
    elseif n==1
        return Polynomial{T}([(a-b)/2, (a+b+2)/2], var)
    end
    p0 = Polynomial{T}([one(T)], var)
    p1 = Polynomial{T}([(a-b)/2, (a+b+2)/2], var)
    px = Polynomial{T}([zero(T), one(T)], var)
    
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
    poly_djacobi(n, [::Type{T}=Int, [var=:x]])

Create a derivative of Jacobi polynomial of order n `Polynomial` object.

Instead of computing the polynomial at a specific point, 
use `Polynomials` package to compute the actual polynomial

Parameters:

 * `n` Order of polynomials
 * `a` Weight α of the Jacobi polynomial
 * `b` Weight β of the Jacobi polynomial
 * `T` Type of polynomial coefficients
 * `var` symbol to be used as variable by the polynomial

# Examples
```julia-repl
julia> p = poly_djacobi(4, 1, 1)
Polynomials.Polynomial(-17.5*x + 52.5*x^3)

julia> p2 = poly_djacobi(4, 1, 1, Rational{Int})
Polynomials.Polynomial(-35//2*x + 105//2*x^3)

julia> p3 = poly_djacobi(4, 1, 1, BigFloat, :y)
Polynomials.Polynomial(-17.5*y + 52.5*y^3)
```
"""
poly_djacobi(n, a, b, ::Type{T}=Float64, var=:x) where {T<:Number} = derivative(poly_jacobi(n,a,b,T,var))
