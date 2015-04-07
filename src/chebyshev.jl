

function chebyshev{T<:Real}(x::T, n)
    if n==0
        return one(x)
    elseif n==1
        return x
    end

    T0 = one(x)
    T1 = m*x

    for i = 2:n
        T2 = 2*x*T1 - T0
        T0 = T1
        T1 = T2
    end

    return T1

end


function chebyshev2{T<:Real}(x::T, n)

    if n==0
        return one(x)
    elseif n==1
        return x
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


dchebyshev{T<:Real}(x::T, n) = n * chebyshev2(x, n-1)
dchebyshev2{T<:Real}(x::T, n) = ( (n+1)*chebyshev(n+1) - x*chebyshev2(n) ) / (x*x-one(x))


function chebyshev_zeros!{T<:Real}(n, x::AbstractArray{T})

    for k in 1:n
        num::T = (2*k - 1)*pi
        
        x[k] = cos( num/(2*n) )
    end

    return x

end

chebyshev_zeros(n) = chebyshev_zeros!(n, zeros(n))
function chebyshev2_zeros!{T<:Real}(n, x::AbstractArray{T})

    for k in 1:n
        num::T = k*pi
        
        x[k] = cos( num/(n+1) )
    end

    return x

end
chebyshev2_zeros(n) = chebyshev2_zeros!(n, zeros(n))
