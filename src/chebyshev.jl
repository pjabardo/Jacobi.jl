
"""
Compute the `n-th` degree Chebyshev polynomial of first and second kind and its derivatives

 * `chebyshev(x, n)` Chebyshev polynomial of the first kind
 * `chebyshev2(x, n)` Chebyshev polynomial of the second kind
 * `dchebyshev(x, n)` Derivative of Chebyshev polynomial of the first kind
 * `dchebyshev2(x, n)` Derivative of Chebyshev polynomial of the second kind

There are also functions that compute the zeros of Chebyshev polynomials:

 * `chebyshev_zeros`
 * `chebyshev2_zeros`
 * `chebyshev_zeros!`
 * `chebyshev2_zeros!`

# Example
```julia-repl

julia> chebyshev(0.3, 5)
0.9988800000000001

julia> chebyshev2(0.3, 5)
1.01376

julia> dchebyshev(0.3, 5)
0.24800000000000044

julia> dchebyshev2(0.3, 5)
-1.3439999999999992

julia> chebyshev_zeros(5)
5-element Array{Float64,1}:
 -0.9510565162951535
 -0.5877852522924731
 -0.0
  0.587785252292473
  0.9510565162951536

julia> chebyshev_zeros(5, BigFloat)
5-element Array{BigFloat,1}:
 -0.9510565162951535721164393333793821434056986341257502224473056444301531700851959
 -0.5877852522924731291687059546390727685976524376431459910722724807572784741623414
 -0.0
  0.5877852522924731291687059546390727685976524376431459910722724807572784741623414
  0.9510565162951535721164393333793821434056986341257502224473056444301531700851873


"""
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


dchebyshev(x, n) = n * chebyshev2(x, n-1)



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




function chebyshev_zeros!(n, x::AbstractArray{T}) where {T<:Number}

    for k in 1:n
        num::T =  (2*k - 1)
        
        x[k] = -cospi( num/(2*n) )
    end

    return x

end

chebyshev_zeros(n, ::Type{T}=Float64) where {T<:Number} = chebyshev_zeros!(n, zeros(T, n))

function chebyshev2_zeros!(n, x::AbstractArray{T}) where {T<:Number}

    for k in 1:n
        num::T = k
        
        x[k] = -cospi( num/(n+1) )
    end

    return x

end
chebyshev2_zeros(n, ::Type{T}=Float64) where {T<:Number} = chebyshev2_zeros!(n, zeros(T, n))

@doc (@doc chebyshev) chebyshev2

@doc (@doc chebyshev) dchebyshev
@doc (@doc chebyshev) dchebyshev2

@doc (@doc chebyshev) chebyshev_zeros!
@doc (@doc chebyshev) chebyshev2_zeros!
@doc (@doc chebyshev) chebyshev_zeros
@doc (@doc chebyshev) chebyshev2_zeros

