using Polynomials

function poly_chebyshev{T<:Number}(n, ::Type{T}=Int, var=:x)

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

poly_dchebyshev{T<:Number}(n, ::Type{T}=Int, var=:x) = polyder(chebyshev(n,T, var))



function poly_chebyshev2{T<:Number}(n, ::Type{T}=Int, var=:x)

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

poly_dchebyshev2{T<:Number}(n, ::Type{T}=Int, var=:x) = polyder(chebyshev2(n,T, var))


function poly_legendre{T<:Number}(n, ::Type{T}=Float64, var=:x)


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


poly_dlegendre{T<:Number}(n, ::Type{T}=Int, var=:x) = polyder(legendre(n,T, var))





function poly_jacobi{T<:Number}(n, a, b, ::Type{T}=Float64, var=:x)

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


poly_djacobi{T<:Number}(n, a, b, ::Type{T}=Float64, var=:x) = polyder(jacobi(n,a,b,T,var))
