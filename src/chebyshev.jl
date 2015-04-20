

function chebyshev(x, n)
    if n==0
        return one(x)
    elseif n==1
        return x
    end

    T0 = one(x)
    T1 = x

    for i = 2:n
        T2 = 2*x*T1 - T0
        T0 = T1
        T1 = T2
    end

    return T1

end

function chebyshev!{T<:Real}(x::AbstractArray{T}, n, y::AbstractArray{T})

    m = length(x)
    for i = 1:m
        y[i] = chebyshev(x[i], n)
    end
    return y
end

chebyshev{T<:Real}(x::AbstractArray{T}, n) = chebyshev!(x, n, zeros(x))

function chebyshev2(x, n)

    if n==0
        return one(x)
    elseif n==1
        return 2x
    end

    U0 = one(x)
    U1 = 2*x

    for i = 2:n
        U2 = 2*x*U1 - U0
        U0 = U1
        U1 = U2
    end

    return U1

    
end

function chebyshev2!{T<:Real}(x::AbstractArray{T}, n, y::AbstractArray{T})

    m = length(x)
    for i = 1:m
        y[i] = chebyshev2(x[i], n)
    end
    return y
end

chebyshev2{T<:Real}(x::AbstractArray{T}, n) = chebyshev2!(x, n, zeros(x))


dchebyshev(x, n) = n * chebyshev2(x, n-1)
function dchebyshev!{T<:Real}(x::AbstractArray{T}, n, y::AbstractArray{T})

    m = length(x)
    for i = 1:m
        y[i] = dchebyshev(x[i], n)
    end
    return y
end

dchebyshev{T<:Real}(x::AbstractArray{T}, n) = dchebyshev!(x, n, zeros(x))


function dchebyshev2(x, n)

    if n == 0
        return zero(x)
    elseif n==1
        return 2*one(x)
    end

    dU0 = zero(x)
    dU1 = 2
    U0 = one(x)
    U1 = 2*x

    for i = 2:n
        U2 = 2*x*U1 - U0
        dU2 = 2*U1 + 2*x*dU1 - dU0
        
        U0 = U1
        U1 = U2
        dU0 = dU1
        dU1 = dU2
    end

    return dU1
end




function dchebyshev2!{T<:Real}(x::AbstractArray{T}, n, y::AbstractArray{T})

    m = length(x)
    for i = 1:m
        y[i] = dchebyshev2(x[i], n)
    end
    return y
end

dchebyshev2{T<:Real}(x::AbstractArray{T}, n) = dchebyshev2!(x, n, zeros(x))




function chebyshev_zeros!{T<:FloatingPoint}(n, x::AbstractArray{T})

    for k in 1:n
        num::T =  (2*k - 1)*pi
        
        x[k] = -cos( num/(2*n) )
    end

    return x

end

chebyshev_zeros{T<:FloatingPoint}(n, ::Type{T}=Float64) = chebyshev_zeros!(n, zeros(T, n))

function chebyshev2_zeros!{T<:FloatingPoint}(n, x::AbstractArray{T})

    for k in 1:n
        num::T = k*pi
        
        x[k] = -cos( num/(n+1) )
    end

    return x

end
chebyshev2_zeros{T<:FloatingPoint}(n, ::Type{T}=Float64) = chebyshev2_zeros!(n, zeros(T, n))

