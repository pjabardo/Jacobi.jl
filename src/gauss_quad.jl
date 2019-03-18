using SpecialFunctions

abstract type QUADRATURE_TYPE end

struct GJ <: QUADRATURE_TYPE end
struct GLJ <: QUADRATURE_TYPE end
struct GRJM <: QUADRATURE_TYPE end
struct GRJP <: QUADRATURE_TYPE end

"""
Gauss-type quadrature

Numerical integrals in the domain [-1,1] can be computed
from knowledge of the function in a set of nodes and the
corresponding nodes, such that

```math
\\int_{-1}^1 (1-x)^a (1+x)^b (1-xf(x)\\:dx \\approx \\sum_{i=1}^N w^{a,b}_i f(x_i)
```

The parameters `a` and `b` form a famility of quadrature rules. But if one or both
of the ends are specified, other quadrature families are possible:

 * No ends are specified, Gauss-Jacobi quadrature (GJ)
 * Both ends are specified, Gauss-Lobatto-Jacobi quadrature (GLJ)
 * A single end is specified, Gauss-Radau-Jacobi quadrature (if +1 is specified, GRJP, if -1, GRJM)

To compute the nodes, the following functions are available:

 * `zgj` (Gauss-Jacobi)
 * `zglj` (Gauss-Lobatto-Jacobi)
 * `zgrjm` (Gauss-Radau-Jacobi, -1)
 * `zgrjp` (Gauss-Radau-Jacobi, +1)

All these functions have the following parameters:

 * `Q`  Number of quadrature nodes
 * `a` (a) weight
 * `b` (b) weight
 *  type Data type to be used, `Float64` is the default

To compute the weights, first the zeros (`z`) should be computed and then the weights are computed
with the following functions:

 * `wgj(z, a, b)`
 * `wglj(z, a, b)`
 * `wgrjm(z, a, b)`
 * `wgrjp(z, a, b)`

 ### Derivatives

The nodes used in the quadrature rules are convenient when using high order Lagrange interpolation
To compute derivatives. The following functions are used to compute the derivative matrix such that
`du = D*u`:

 * `dgj(z, a, b)`
 * `dglj(z, a, b)`
 * `dgrjm(z, a, b)`
 * `dgrjp(z, a, b)`

### Examples

See the notebooks availbale with the package.

"""
function zgj(Q, a, b, ::Type{T}=Float64) where {T<:Number}
    jacobi_zeros(Q, a, b, T)
end

zgj(Q) = zgj(Q, 0.0, 0.0)
zgj(Q, a) = zgj(Q, a, zero(a))

function zglj(Q, a, b, ::Type{T}=Float64) where {T<:Number}
    z = jacobi_zeros(Q-2, a+1, b+1, T)
    o = one(T)
    return [-o; z; o]
end
zglj(Q) = zglj(Q, 0.0, 0.0)
zglj(Q, a) = zglj(Q, a, zero(a))


function zgrjm(Q, a, b, ::Type{T}=Float64) where {T<:Number}
    z = jacobi_zeros(Q-1, a, b+1, T)
    return [-one(T); z]
end
zgrjm(Q) = zgrjm(Q, 0.0, 0.0)
zgrjm(Q, a) = zgrjm(Q, a, zero(a))


function zgrjp(Q, a, b, ::Type{T}=Float64) where {T<:Number}
    z = jacobi_zeros(Q-1, a+1, b, T)
    return [z; one(T)]
end
zgrjp(Q) = zgrjp(Q, 0.0, 0.0)
zgrjp(Q, a) = zgrjp(Q, a, zero(a))


function wgj(z::AbstractArray{T}, alpha=0, beta=0) where {T<:Number}
    a = convert(T, alpha)
    b = convert(T, beta)

    Q::Int = length(z)
    o = one(T)
    coef = 2^(a+b+1) * gamma_ratio(a+Q+1, Q+o) * gamma_ratio(b+Q+1, a+b+Q+1)
    w = [djacobi(zz, Q, a, b) for zz=z]

    for i = 1:Q
        ww = w[i]
        x = z[i]
        w[i] = o / (ww*ww) * coef / (o - x*x)
    end

    return w
end


function wglj(z::AbstractArray{T}, alpha=0, beta=0) where {T<:Number}
    a = convert(T, alpha)
    b = convert(T, beta)
    o = one(T)
    Q = length(z)

    coef = 2^(a+b+1) / (Q-o) * gamma_ratio(a+Q, Q*o) * gamma_ratio(b+Q, a+b+Q+1)
    
    w = [jacobi(zz, Q-1, a, b) for zz=z]
    w[1] = (b+1) * coef / (w[1]*w[1])
    w[Q] = (a+1) * coef / (w[Q]*w[Q])

    for i = 2:(Q-1)
        ww = w[i]
        w[i] = coef / (ww * ww)
    end

    return w
end

function wgrjm(z::AbstractArray{T}, alpha=0, beta=0) where {T<:Number}
    a = convert(T, alpha)
    b = convert(T, beta)
    o = one(T)

    Q = length(z)

    coef = 2^(a+b) / (b+Q) * gamma_ratio(a+Q, Q*o) * gamma_ratio(b+Q, a+b+Q+1)

    w = [jacobi(zz, Q-1, a, b) for zz=z]

    for i = 1:Q
        ww = w[i]
        w[i] = coef / (ww*ww) * (o - z[i])
    end

    w[1] *= (b + o)

    return w
end


function wgrjp(z::AbstractArray{T,1}, alpha=0, beta=0) where {T<:Number}
    a = convert(T, alpha)
    b = convert(T, beta)
    Q = length(z)
    o = one(T)

    coef = 2^(a+b) / (a+Q) * gamma_ratio(a+Q, Q*o) * gamma_ratio(b+Q, a+b+Q+1)

    w = [jacobi(zz, Q-1, a, b) for zz=z]

    for i = 1:Q
        ww = w[i]
        w[i] = coef / (ww*ww) * (o + z[i])
    end

    w[Q] *= (a + o)

    return w
end


function dgj(z::AbstractArray{T,1}, alpha=0, beta=0) where {T<:Number}

    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    o = one(T)

    djac = [djacobi(zz, Q, alpha, beta) for zz=z]

    D = zeros(T, Q, Q)
    for i = 1:Q
        for k = 1:Q
            if i != k
                D[i,k] = (djac[i]/djac[k]) / (z[i]-z[k])
            else
                D[i,i] = (a-b + (a + b + 2) * z[i]) / (o - z[k]^2) / 2
            end
        end
    end

    return D
end


function dglj(z::AbstractArray{T,1}, alpha=0, beta=0) where {T<:Number}

    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    o = one(T)

    djac = zeros(T,Q)
    djac[1] = (-1)^Q * 2 * gamma_ratio(Q+b, Q-o) / gamma(b+2)
    djac[Q] = -2 * gamma_ratio(Q+a, Q-o) / gamma(a+2)
    for i = 2:(Q-1)
        djac[i] = (o-z[i]*z[i]) * djacobi(z[i], Q-2, a+1, b+1)
    end

    D = zeros(T, Q, Q)
    for i = 1:Q
        for k = 1:Q
            if i != k
                D[i,k] = (djac[i]/djac[k]) / (z[i]-z[k])
            else
                D[i,i] = (a - b + (a + b)*z[i]) / (2*(o - z[i]^2))
            end
        end
    end

    D[1,1] = (a - (Q-1)*(Q+a+b)) / (2*(b + 2))
    D[Q,Q] = -(b - (Q-1)*(Q+a+b)) / (2*(a + 2))

    return D

end


function dgrjm(z::AbstractArray{T,1}, alpha=0, beta=0) where {T<:Number}

    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    o = one(T)

    djac = zeros(T,Q)
    for i = 2:Q
        djac[i] = (1+z[i]) * djacobi(z[i], Q-1, a, b+1)
    end
    djac[1] = (-1)^(Q-1) * gamma_ratio(Q+b+1, Q*o) / gamma(b+2)

    D = zeros(T, Q, Q)
    for i = 1:Q
        for k = 1:Q
            if i != k
                D[i,k] = (djac[i]/djac[k]) / (z[i]-z[k])
            else
                D[i,i] = (a - b + 1 + (a + b + 1)*z[i]) / (2*(1 - z[i]^2))
            end
        end
    end

    D[1,1] = -(Q-1)*(Q+a+b+1) / (2*(b + 2))

    return D

end


function dgrjp(z::AbstractArray{T,1}, alpha=0, beta=0) where {T<:Number}

    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    o = one(T)

    djac = zeros(T,Q)
    for i = 1:(Q-1)
        djac[i] = (1-z[i]) * djacobi(z[i], Q-1, a+1, b)
    end
    djac[Q] = - gamma_ratio(Q+a+1, Q*o) / gamma(a+2)

    D = zeros(T, Q, Q)
    for i = 1:Q
        for k = 1:Q
            if i != k
                D[i,k] = (djac[i]/djac[k]) / (z[i]-z[k])
            else
                D[i,i] = (a - b - 1 + (a + b + 1)*z[i]) / (2*(1 - z[i]^2))
            end
        end
    end

    D[Q,Q] = (Q-1)*(Q+a+b+1) / (2*(a + 2))

    return D

end


"""
Abstract interface of Gauss-type quadrature rules

Can be used for any `AbstractFloat` type data.
"""
struct Quadrature{T<:Number,QT<:QUADRATURE_TYPE}
    "Number of quadrature nodes"
    Q::Int
    "a weight"
    a::T
    "b weight"
    b::T
    "Quadrature nodes"
    z::Array{T, 1}
    "Quadrature weights"
    w::Array{T, 1}
    "Quadrature derivative matrix"
    D::Array{T, 2}
end

"""
Return the zeros of a Gauss type quadrature
"""
qzeros(::Type{GJ}, Q, a=0, b=0, ::Type{T}=Float64) where {T<:Number} = zgj(Q, a, b, T)
qzeros(::Type{GLJ}, Q, a=0, b=0, ::Type{T}=Float64) where {T<:Number} = zglj(Q, a, b, T)
qzeros(::Type{GRJM}, Q, a=0, b=0, ::Type{T}=Float64) where {T<:Number} = zgrjm(Q, a, b, T)
qzeros(::Type{GRJP}, Q, a=0, b=0, ::Type{T}=Float64) where {T<:Number} = zgrjp(Q, a, b, T)

"""
Return the weights of a Gauss type quadrature
"""
qweights(::Type{GJ}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = wgj(z, a, b)
qweights(::Type{GLJ}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = wglj(z, a, b)
qweights(::Type{GRJM}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = wgrjm(z, a, b)
qweights(::Type{GRJP}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = wgrjp(z, a, b)

"""
Return the derivative matrix of a Gauss type quadrature
"""
qdiff(::Type{GJ}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = dgj(z, a, b)
qdiff(::Type{GLJ}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = dglj(z, a, b)
qdiff(::Type{GRJM}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = dgrjm(z, a, b)
qdiff(::Type{GRJP}, z::AbstractArray{T}, a=0, b=0) where {T<:Number} = dgrjp(z, a, b)


"""
Create a `Quadrature` object given its type, order and weights.
"""
function Quadrature(::Type{QT}, Q, a=0, b=0,
                    ::Type{T}=Float64) where {T<:Number, QT<:QUADRATURE_TYPE}
    aa = convert(T, a)
    bb = convert(T, b)
    z = qzeros(QT, Q, aa, bb, T)
    w = qweights(QT, z, aa, bb)
    D = qdiff(QT, z, aa, bb)
    Quadrature{T,QT}(Q, aa, bb, z, w, D)
end

"Return quadrature type"
qtype(q::Quadrature{T,QT}) where {T,QT} = QT
"Return quadrature nodes"
qzeros(q::QT) where {QT<:Quadrature} = q.z
"Return quadrature weights"
qweights(q::QT) where {QT<:Quadrature} = q.w
"Return quadrature derivative matrix"
qdiff(q::QT) where {QT<:Quadrature} = q.D
"Return number of quadrature nodes"
num_points(q::QT) where {QT<:Quadrature} = q.Q
"Return quadrature `a` weight"
qalpha(q::QT) where {QT<:Quadrature} = q.a
"Return quadrature `b` weight"
qbeta(q::QT) where {QT<:Quadrature} = q.b

"""
Compute the Lagrange polynomial

This function computes the Lagrange polynomial at point `x` corresponding to
the `i`-th node of the set `z`

There is also a modifying version `lagrange!` used to computing the polynomials
at several points of an array `x`
"""
function lagrange(i, x, z)
    nz = length(z)

    l = one(z[1])

    for k = 1:(i-1)
        l = l * (x-z[k]) / (z[i]-z[k])
    end

    for k = (i+1):nz
        l = l * (x-z[k]) / (z[i]-z[k])
    end

    return l
end



"""
Interpolation matrix

Often it is necessary to compute the values of a function approximated using
Lagrange interpolation through a set of points `z`. If this will be repeated
often, a matrix can be computed that allows the easy computation using the simple
expression `fx = Imat * fz`
"""
function interp_mat(x::AbstractArray{T}, z::AbstractArray{T}) where {T<:Number}

    Q = length(z)
    np = length(x)

    Imat = zeros(T, np, Q)
    for i = 1:Q, k=1:np
        Imat[k,i] = lagrange(i, x[k], z)
    end

    return(Imat)

end




@doc (@doc zgj) zglj
@doc (@doc zgj) zgrjm
@doc (@doc zgj) zgrjp

@doc (@doc zgj) wgj
@doc (@doc zgj) wglj
@doc (@doc zgj) wgrjm
@doc (@doc zgj) wgrjp

@doc (@doc zgj) dgj
@doc (@doc zgj) dglj
@doc (@doc zgj) dgrjm
@doc (@doc zgj) dgrjp

