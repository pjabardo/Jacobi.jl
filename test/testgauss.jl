### Testing quadrature related function

Nmax= 20
a = rand(1:6, Nmax)  # Coefficient for a ninth degree polynomial

N = 10
Q = 5
p = Poly(a[1:N])
p_int = polyint(p)
X = polyval(p_int, 1.0) - polyval(p_int, -1.0)
q = Quadrature(GJ, Q, 0.0, 0.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)

# I'm not getting the precision I need for BigFloats. Should investigate further
#X = polyval(p1_int, BigFloat(1.0)) - polyval(p1_int, BigInt(-1.0))
#q2 = Quadrature(GJ, 5, BigFloat(0), BigFloat(0), BigFloat)
#y = polyval(p1, qzeros(q2))
#Xg= sum(qweights(q2) .* y)
#@test abs(Xg-X) 0.0 200*eps(Xg)


# Weights a=1, b=2
p2 = p * Poly([1, -1]) * Poly([1, 1])^2
p2_int = polyint(p2)
X = polyval(p2_int, 1.0) - polyval(p2_int, -1.0)
q = Quadrature(GJ, Q, 1.0, 2.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)

# Test Gauss-Radau:
# GRJM
N = 9
Q = 5
p = Poly(a[1:N])
p_int = polyint(p)
X = polyval(p_int, 1.0) - polyval(p_int, -1.0)
q = Quadrature(GRJM, Q, 0.0, 0.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)

# GRJP
q = Quadrature(GRJP, Q, 0.0, 0.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)


# Weights a=1, b=2 GRJM
p2 = p * Poly([1, -1]) * Poly([1, 1])^2
p2_int = polyint(p2)
X = polyval(p2_int, 1.0) - polyval(p2_int, -1.0)
q = Quadrature(GRJM, Q, 1.0, 2.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)


# Weights a=1, b=2 GRJP
p2 = p * Poly([1, -1]) * Poly([1, 1])^2
p2_int = polyint(p2)
X = polyval(p2_int, 1.0) - polyval(p2_int, -1.0)
q = Quadrature(GRJP, Q, 1.0, 2.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)





# Test Gauss-Lobatto
N = 8
Q = 5
p = Poly(a[1:N])
p_int = polyint(p)
X = polyval(p_int, 1.0) - polyval(p_int, -1.0)
q = Quadrature(GLJ, Q, 0.0, 0.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)

# Weights a=1, b=2 GLJ
p2 = p * Poly([1, -1]) * Poly([1, 1])^2
p2_int = polyint(p2)
X = polyval(p2_int, 1.0) - polyval(p2_int, -1.0)
q = Quadrature(GRJM, Q, 1.0, 2.0, Float64)
y = polyval(p, qzeros(q))
Xg = sum(qweights(q) .* y)
@test abs(Xg-X) ≈ 0.0 atol=200*eps(Xg)


# Testing derivatives:
# GJ
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GJ, Q, 0, 0)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))


# GRJM
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GRJM, Q, 0, 0)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))



# GRJP
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GRJP, Q, 0, 0)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))


# GLJ
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GLJ, Q, 0, 0)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))



# Testing with weights a and b
a1 = 0.5
b1 = 0.9

# GJ
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GJ, Q, a1, b1)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))


# GRJM
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GRJM, Q, a1, b1)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))



# GRJP
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GRJP, Q, a1, b1)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))


# GLJ
Q = 10
N = 10
p = Poly(a[1:N])
dp = polyder(p)
q = Quadrature(GLJ, Q, a1, b1)
z = qzeros(q)
y = polyval(p, z)
dy = polyval(dp, z)
dyg = qdiff(q) * y
@test maximum(abs,dy-dyg) ≈ 0.0 atol=1000*eps(maximum(abs,dy))



# Test interpolation:
# GJ

N = 10
Q = 10
p = Poly(a[1:N])
x = Compat.range(-1, stop=1, length=101)
ue = polyval(p, x)
q = Quadrature(GLJ, Q, 1.0, 0.5)
z = qzeros(q)
u = polyval(p, z)
I = interp_mat(x, z)
ui = I * u
@test maximum(abs,ui-ue) ≈ 0.0 atol=200*eps(maximum(abs,u))

