using Base.Test


import Jacobi


z = [0.05, 0.25, 0.50, 0.75, 0.95]

L = [-0.49625  -0.40625  -0.125  0.34375  0.85375;
     -0.0746875  -0.3359375  -0.4375000  -0.0703125  0.7184375;
     0.11857899  0.17682442  -0.26789856  0.31033185  -0.26842182;
     -0.2129835  0.2212002  -0.1882286  0.2643745  -0.3548803]


m = [2, 3, 9, 10]
np  = length(z)

# Testing Legendre polynomials - source Abramowitz
for i = 1:4
    for k = 1:np
        @test L[i,k] ≈ Jacobi.legendre(z[k], m[i]) atol=1e-7
    end
end

# Testing Legendre polynomials comparing to Jacobi polynomials
for i = 1:4
    for k = 1:np
        @test Jacobi.legendre(z[k], m[i]) ≈ Jacobi.jacobi(z[k], m[i], 0, 0) atol=10*eps(1.0)
    end
end


m2 = [3, 4, 9, 10]

dL = [-1.48125 -1.03125 0.3750 2.71875 5.26875;
      -0.3728125 -1.6015625 -1.56250 1.7578125 7.8790625;
      2.1946113 -1.8879372 0.7237244 -0.7228714 10.2440570;
      1.2955205 1.2962599 -2.3171234 2.5611649 7.0476358]

# Testing derivatives of Legendre polynomials:
for i = 1:4
    for k = 1:np
        @test dL[i,k] ≈ Jacobi.dlegendre(z[k], m2[i]) atol=1e-7
    end
end


# Test other derivatives:

calc_deriv(f, x, dx) = (f(x+dx) - f(x-dx)) / (2*dx)


a = 0.8
b = 0.3

for i = 1:4
    for k = 1:np
        d = calc_deriv(x->Jacobi.jacobi(x, m2[i], a, b), z[k], 1e-8)
        @test d ≈ Jacobi.djacobi(z[k], m2[i], a, b) atol=5e-7
    end
end



# Testing Calculation of zeros
      
a = [-1.0, -0.5, 0, 0.5, 1.0]
b = [-0.5, 0.5, 0, 0.5, -0.5]

na = length(a)

m = [3:5;]

nm = length(m)

for i = 1:na
    for k = 1:nm
        Q = m[k]
        z = Jacobi.jacobi_zeros(Q, a[i], b[i])
        for j = 1:Q
            y = Jacobi.jacobi(z[j], Q, a[i], b[i])
            @test abs(y) < 500*eps(1.0)
        end
    end
end

# Testing for bigfloats:
for i = 1:na
    for k = 1:nm
        Q = m[k]
        zb = Jacobi.jacobi_zeros(Q, a[i], b[i], BigFloat)
        for j = 1:Q
            yb = Jacobi.jacobi(zb[j], Q, a[i], b[i])
            @test yb ≈ BigFloat(0) atol=500*eps(BigFloat(1))
        end
    end
end
            
            


function peval(x, p)
    n = length(p)
    y = one(x)*p[n]
    for i = (n-1):-1:1
        y = p[i] + x*y
    end
    return y
end
function peval(x::AbstractArray, p)
    y = zeros(x)
    for i = 1:length(x)
        y[i] = peval(x[i], p)
    end
    return y
end

function dpoly(p)
    n = length(p)
    if (n==1)
        return([0*p[1]])
    end
    a = zeros(eltype(p), n-1)
    for i = 1:(n-1)
        a[i] = p[i+1]*i
    end
    return a
end

function ipoly(p)
    n = length(p)
    a = zeros(eltype(p), n+1)
    a[1] = zero(eltype(p))
    for i = 1:n
        a[i+1] = a[i] / n
    end
    return a
end


  

       
# Chebyshev polynomials
# Coefficients from Abramowiz

xx = linspace(-1.0, 1.0, 21)

using Polynomials

t11 = Poly([0, -11, 0, 220, 0, -1232, 0, 2816, 0, -2816, 0, 1024])
u11 = Poly([0, -12, 0, 280, 0, -1792, 0, 4608, 0, -5120, 0, 2048])

pt11 = Jacobi.poly_chebyshev(11)
pu11 = Jacobi.poly_chebyshev2(11)
pdt11 = Jacobi.poly_dchebyshev(11)
pdu11 = Jacobi.poly_dchebyshev2(11)

dt11 = polyder(t11)
du11 = polyder(u11)

@test t11 == pt11
@test u11 == pu11
@test dt11 == pdt11
@test du11 == pdu11


y1 = zeros(xx)

# Testing Chebyshev polynomials of the first kind
Jacobi.chebyshev!(xx, 11, y1)
y2 = polyval(t11, xx)
@test maximum(abs, y1-y2) < 1e-12

# Testing Chebyshev polynomials of the second kind
Jacobi.chebyshev2!(xx, 11, y1)
y2 = polyval(u11, xx)
@test maximum(abs, y1-y2) < 1e-12


# Testing Chebyshev polynomials of the first kind
Jacobi.dchebyshev!(xx, 11, y1)
y2 = polyval(dt11, xx)
@test maximum(abs, y1-y2) < 1e-12

# Testing Chebyshev polynomials of the second kind
Jacobi.dchebyshev2!(xx, 11, y1)
y2 = polyval(du11, xx)
@test maximum(abs, y1-y2) < 3e-12


# Testing legendre polynomials:
leg11 = Poly([0, -693, 0, 15015, 0, -90090, 0, 218790, 0, -230945, 0, 88179])/256
pleg11 = Jacobi.poly_legendre(11)
@test maximum(abs, leg11.a[1:12] - pleg11.a[1:12]) < 200*eps(500.0)



# Testing Chebyshev zeros:
for k = 1:20
    z  = Jacobi.chebyshev_zeros(k)
    z1 = Jacobi.jacobi_zeros(k, -0.5, -0.5)
    @test maximum(abs, z-z1) < 1e-13
end


# Testing array function
y1 = zeros(xx)
y2 = zeros(xx)


Jacobi.jacobi!(xx, 5, 0.8, 0.3, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.jacobi(xx[i], 5, 0.8, 0.3)
end
@test maximum(abs, y1-y2) < eps(Float64)


Jacobi.djacobi!(xx, 5, 0.8, 0.3, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.djacobi(xx[i], 5, 0.8, 0.3)
end
@test maximum(abs, y1-y2) < eps(Float64)



Jacobi.legendre!(xx, 5, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.legendre(xx[i], 5)
end
@test maximum(abs, y1-y2) < eps(Float64)


Jacobi.dlegendre!(xx, 5, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.dlegendre(xx[i], 5)
end
@test maximum(abs, y1-y2) < eps(Float64)


Jacobi.chebyshev!(xx, 5, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.chebyshev(xx[i], 5)
end
@test maximum(abs, y1-y2) < eps(Float64)


Jacobi.chebyshev2!(xx, 5, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.chebyshev2(xx[i], 5)
end
@test maximum(abs, y1-y2) < eps(Float64)

Jacobi.dchebyshev!(xx, 5, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.dchebyshev(xx[i], 5)
end
@test maximum(abs, y1-y2) < eps(Float64)


Jacobi.dchebyshev2!(xx, 5, y1)
for i in 1:length(xx)
    y2[i] = Jacobi.dchebyshev2(xx[i], 5)
end
@test maximum(abs, y1-y2) < eps(Float64)





true
