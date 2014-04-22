

function zgj{T<:Real}(Q, a::T, b::T)
    jacobi_zeros(Q, a, b)
end
zgj(Q) = zgj(Q, 0.0, 0.0)
zgj(Q, a) = zgj(Q, a, zero(a))

function zglj{T<:Real}(Q, a::T, b::T)
    z = jacobi_zeros(Q-2, a+1, b+1)
    o = one(T)
    return [-o, z, o]
end
zglj(Q) = zglj(Q, 0.0, 0.0)
zglj(Q, a) = zglj(Q, a, zero(a))


function zgrjm{T<:Real}(Q, a::T, b::T)
    z = jacobi_zeros(Q-1, a, b+1)
    return [-1, z]
end
zgrjm(Q) = zgrjm(Q, 0.0, 0.0)
zgrjm(Q, a) = zgrjm(Q, a, zero(a))


function zgrjp{T<:Real}(Q, a::T, b::T)
    z = jacobi_zeros(Q-1, a+1, b)
    return [z, one(T)]
end
zgrjp(Q) = zgrjp(Q, 0.0, 0.0)
zgrjp(Q, a) = zgrjp(Q, a, zero(a))


function wgj{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    a = convert(T, alpha)
    b = convert(T, beta)

    Q::Int = length(z)
    
    coef = 2^(a+b+1) * ( gamma(a+Q+1) / factorial(Q) ) * (gamma(b+Q+1) / gamma(a+b+Q+1))

    w = [djacobi(zz, Q, a, b) for zz=z]
    for i = 1:Q
        ww = w[i]
        x = z[i]
        w[i] = 1 / (ww*ww) * coef / (1 - x*x)
    end

    return w
end


function wglj{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    a = convert(T, alpha)
    b = convert(T, beta)

    Q = length(z)

    coef = 2^(a+b+1) / (Q-1) * (gamma(a+Q) / factorial(Q-1)) * (gamma(b+Q) / gamma(a+b+Q+1))
    
    w = [jacobi(zz, Q-1, a, b) for zz=z]
    w[1] = (b+1) * coef / (w[1]*w[1])
    w[Q] = (a+1) * coef / (w[Q]*w[Q])

    for i = 2:(Q-1)
        ww = w[i]
        w[i] = coef / (ww * ww)
    end

    return w
end

function wgrjm{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    a = convert(T, alpha)
    b = convert(T, beta)
    
    Q = length(z)

    coef = 2^(a+b) / (b+Q) * (gamma(a+Q) / factorial(Q-1)) * (gamma(b+Q) / gamma(a+b+Q+1))

    w = [jacobi(zz, Q-1, a, b) for zz=z]

    for i = 1:Q
        ww = w[i]
        w[i] = coef / (ww*ww) * (1 - z[i])
    end

    w[1] *= (b + 1)

    return w
end


function wgrjp{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    a = convert(T, alpha)
    b = convert(T, beta)
    Q = length(z)

    coef = 2^(a+b) / (a+Q) * (gamma(a+Q) / factorial(Q-1)) * (gamma(b+Q) / gamma(a+b+Q+1))

    w = [jacobi(zz, Q-1, a, b) for zz=z]

    for i = 1:Q
        ww = w[i]
        w[i] = coef / (ww*ww) * (1 + z[i])
    end

    w[Q] *= (a + 1)

    return w
end


function dgj{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    
    djac = [djacobi(zz, Q, alpha, beta) for zz=z]
    
    D = zeros(T, Q, Q)
    for i = 1:Q
        for k = 1:Q
            if i != k
                D[i,k] = (djac[i]/djac[k]) / (z[i]-z[k])
            else
                D[i,i] = (a-b + (a + b + 2) * z[i]) / (1 - z[k]^2) / 2
            end
        end
    end
    
    return D
end


function dglj{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)

    djac = zeros(T,Q)
    djac[1] = (-1)^Q * 2 * gamma(Q + b) / (gamma(Q-1) * gamma(b+2))
    djac[Q] = -2*gamma(Q+a) / (gamma(Q-1)*gamma(a+2))
    for i = 2:(Q-1)
        djac[i] = (1-z[i]*z[i]) * djacobi(z[i], Q-2, a+1, b+1)
    end

    D = zeros(T, Q, Q)
    for i = 1:Q
        for k = 1:Q
            if i != k
                D[i,k] = (djac[i]/djac[k]) / (z[i]-z[k])
            else
                D[i,i] = (a - b + (a + b)*z[i]) / (2*(1 - z[i]^2))
            end
        end
    end
 
    D[1,1] = (a - (Q-1)*(Q+a+b)) / (2*(b + 2))
    D[Q,Q] = -(b - (Q-1)*(Q+a+b)) / (2*(a + 2))

    return D

end


function dgrjm{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)

    djac = zeros(T,Q)
    for i = 2:Q
        djac[i] = (1+z[i]) * djacobi(z[i], Q-1, a, b+1)
    end
    djac[1] = (-1)^(Q-1) *  gamma(Q + b + 1) / (gamma(Q) * gamma(b+2))

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


function dgrjp{T<:Real}(z::Array{T,1}, alpha=0, beta=0)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)

    djac = zeros(T,Q)
    for i = 1:(Q-1)
        djac[i] = (1-z[i]) * djacobi(z[i], Q-1, a+1, b)
    end
    djac[Q] = - gamma(Q + a + 1) / (gamma(Q) * gamma(a+2))

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



function lgj{T<:Real}(i::Int, zz::T, z::Array{T,1}, alpha=0, beta=0, epsmult=50)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    err = epsmult * eps(T)

    if abs(z[i]-zz) < err
        return one(zz)
    end

    return jacobi(zz, Q, a, b) / (djacobi(z[i], Q, a, b) * (zz - z[i]))
end


function lglj{T<:Real}(i::Int, zz::T, z::Array{T,1}, alpha=0, beta=0, epsmult=50)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    err = epsmult * eps(T)
    zi = z[i]
    
    if abs(zi-zz) < err
        return one(zz)
    end
     
    local num::T = (1-zz*zz) * jacobi(zz, Q-2, a+1, b+1)
    local den1::T = -2*zi*jacobi(zi, Q-2, a+1, b+1)
    local den2::T = (1-zi*zi) * djacobi(zi, Q-2, a+1, b+1)

    return num / ((den1 + den2)*(zz-zi))
    
end


function lgrjm{T<:Real}(i::Int, zz::T, z::Array{T,1}, alpha=0, beta=0, epsmult=200)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    err = epsmult * eps(T)
    zi = z[i]
    
    if abs(zi-zz) < err
        return one(zz)
    end
     

    return (1+zz)*jacobi(zz, Q-1, a, b+1) / ( (jacobi(zi, Q-1, a, b+1) + (1+zi)*djacobi(zi, Q-1, a, b+1)) * (zz-zi))
    
    
end



function lgrjp{T<:Real}(i::Int, zz::T, z::Array{T,1}, alpha=0, beta=0, epsmult=200)
    
    Q = length(z)
    a = convert(T, alpha)
    b = convert(T, beta)
    err = epsmult * eps(T)
    zi = z[i]
    
    if abs(zi-zz) < err
        return one(T)
    end
     

    return (1-zz)*jacobi(zz, Q-1, a+1, b) / ( (-jacobi(zi, Q-1, a+1, b) + (1-zi)*djacobi(zi, Q-1, a+1, b)) * (zz-zi))
    
    
end





function igj{T<:Real}(zp::Array{T,1}, z::Array{T,1}, alpha=0, beta=0)

    Q = length(z)
    np = length(zp)
    
    Imat = zeros(T, np, Q)
    for i = 1:Q
        for k = 1:np
            Imat[k,i] = lagrange_gj(i, zp[k], z, alpha, beta)
        end
    end

    return(Imat)

end


function iglj{T<:Real}(zp::Array{T,1}, z::Array{T,1}, alpha=0, beta=0)

    Q = length(z)
    np = length(zp)
    
    Imat = zeros(T, np, Q)
    for i = 1:Q
        for k = 1:np
            Imat[k,i] = lagrange_glj(i, zp[k], z, alpha, beta)
        end
    end

    return(Imat)

end

function igrjm{T<:Real}(zp::Array{T,1}, z::Array{T,1}, alpha=0, beta=0)

    Q = length(z)
    np = length(zp)
    
    Imat = zeros(T, np, Q)
    for i = 1:Q
        for k = 1:np
            Imat[k,i] = lagrange_grjm(i, zp[k], z, alpha, beta)
        end
    end

    return(Imat)

end


function igrjp{T<:Real}(zp::Array{T,1}, z::Array{T,1}, alpha=0, beta=0)

    Q = length(z)
    np = length(zp)
    
    Imat = zeros(T, np, Q)
    for i = 1:Q
        for k = 1:np
            Imat[k,i] = lagrange_grjp(i, zp[k], z, alpha, beta)
        end
    end

    return(Imat)

end
