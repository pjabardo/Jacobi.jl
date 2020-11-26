# Jacobi

[![Build Status](https://github.com/pjabardo/Jacobi.jl/workflows/CI/badge.svg)](https://github.com/pjabardo/Jacobi.jl/actions)
[![Build Status](https://travis-ci.com/pjabardo/Jacobi.jl.svg?branch=master)](https://travis-ci.com/pjabardo/Jacobi.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/pjabardo/Jacobi.jl?svg=true)](https://ci.appveyor.com/project/pjabardo/Jacobi-jl)
[![Coverage](https://codecov.io/gh/pjabardo/Jacobi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/pjabardo/Jacobi.jl)
[![Coverage](https://coveralls.io/repos/github/pjabardo/Jacobi.jl/badge.svg?branch=master)](https://coveralls.io/github/pjabardo/Jacobi.jl?branch=master)


A library that implements Jacobi polynomials and Gauss 
quadrature related operations.

## Notebook

In the notebooks directory there are tow notebooks, as follows, that shows the basic usage of the library.

 * [Calculation of some orthogonal polynomials and their zeros](http://nbviewer.ipython.org/github/pjabardo/Jacobi.jl/blob/master/notebooks/polynomials.ipynb)
 * [Gauss-type qudrature, interpolation and differentiation](http://nbviewer.ipython.org/github/pjabardo/Jacobi.jl/blob/master/notebooks/gauss-quad.ipynb)

   
## Jacobi polynomials

Jacobi polynomials  $P_n^{a,b}(x)$ are implemented in function 
`jacobi(x, n, a, b)` where `x` is where you want to calculate the 
polynomial of degree `n` with weights `a` and `b`.

To calculate the derivative of jacobi polynomials, use function 
`djacobi(x, n, a, b)`. 

When calculating the weights of Gaussian qudrature, it is necessary to determine
the zeros of Jacobi polynomials. The function `jacobi_zeros(m, a, b)` calculates the 
zeros of a Jacobi polynomial of degree `m` with weights `a` and `b`. Also present 
is the non allocating function `jacobi_zeros!(m, a, b, x)`.

Legendre polynomials are a special case when a and b are zero. 
Function `legendre(x, n)` implements this simpler recurrence relation.

Also available are Chebyshev polynomials and the respective derivatives and zeros:
 * `chebyshev(x, n)` Chebyshev polynomial of the first kind $T_n(x)$.
 * `chebyshev2(x, n)` Chebyshev polynomial of the second kind $U_n(x)$.
 * `chebyshev_zeros(n)` and `chebyshev2_zeros(n)` Roots of Chebyshev polynomials
 * `chebyshev_zeros!(n, x)` and `chebyshev2_zeros!(n, x)` Roots of Chebyshev polynomials

     

## Gauss Quadrature

Gaussian quadrature is commonly used to evaluate integrals numerically. To evaluate 
integrals using Gaussian quadrature, a set of nodes and its corresponding weights must 
be known. In the most basic algorithm, no node can be specified by the user. It is 
often convenient to include one or both ends of the domain and this provides
different quadrature rules. 

The different quadrature rules implemented are:

 * Gauss-Jacobi, where no nodes are specified (extension `gj`).
 * Gauss-Lobatto-Jacobi, where both ends of the domain are included (extension `glj`).
 * Gauss-Radau-Jacobi, where the node -1 (leftmost end of the domain) is included 
   (extension `grjm`).
 * Gauss-Radau-Jacobi, where the node +1 (rightmost end of the domain) is included 
   (extension `grjp`).

These formulations implement the numerical integration of
$$
\int_{-1}^1 (1-x)^a(1+x)^b f(x)dx
$$

The functions `zgj`, `zglj`, `zgrjm` and `zgrjp` calculate the nodes of the quadrature
rule. The corresponding weights are calculated using functions `wgj`, `wglj`,
`wgrjm` and `wgrjp`. Of course in both cases, the weights `a` and `b`must be specified. 

As an example

```
fun(x) = sin(x)
z = zglj(5, 0.0, 0.0)
w = wglj(z, 0.0, 0.0)
f = fun(z)
Ix = sum(w .* f)
```
where `Ix` is the estimate of the integral.

Gaussian quadrature is useful because it allows the exact integration of polynomials
using few nodes (the exact order depends on the quadrature rule cited above).

But there is another advantage: the nodes of the quadrature rule is very convenient for
interpolating functions using high order polynomials and is commonly used in high order
finite element procedures such as hp-FEM or spectral element method. In these
 applications, it is necessary to compute derivatives and interpolate data from 
different grids. 

The functions `dgj`, `dglj`, `dgrjm` and `drjp` calculates the derivative matrix as 
shown in the example that follows:

```
fun(x) = sin(x)
z = zglj(5, 0.0, 0.0)
D = dglj(z, 0.0, 0.0)
f = fun(z)
df = D * f
```
where `df` is an estimate of the derivative of the function at the quadrature nodes.

Another important operation is interpolation. If a function is known at some nodes, 
in this case the quadrature nodes, how can we accurately interpolate the function on
other nodes? Since we know the nodes, Lagrangian interpolation is the best way.
The function `lagrange` implements the standard  definition of the Langrangian 
interpolation.  The example below plots the Lagrangian interpolators of
the Gauss-Lobatto-Jacobi quadrature points for 5 nodes.

```
using PyPlot
Q = 5
z = zglj(Q)
nx = 201
x = -1.0:0.01:1.0
y = zeros(nx, Q)
for k = 1:Q, i=1:nx
  y[i,k] = lagrange(k, x[i], z)
end

for k=1:Q
  plot(x, y[:,k])
end
```

If the operation above is to be repeated often, pre-calculating the Lagrangian 
interpolators is useful and an Interpolation matrix can be calculated. The 
following example illustrates the use of the interpolation matrix that can 
be computed with the function `interp_mat`. 

```
using PyPlot
Q = 5
z = zglj(Q)
nx = 201
x = -1.0:0.01:1.0
ye = sin(pi*z)
ye2 = sin(pi*x)
Im = interp_mat(x, z)
y = Im * ye
plot(z, ye, "o")
plot(x, ye2, "r-")
plot(x, y, "b-")
```
increasing the number of quadrature points the interpolated function (blue line) 
becomes more accurate.


## References

This package was implemented using both references below. 

 * Spectral/hp Element Methods for CFD, 2nd edition, Karniadakis and Sherwin, 2005.
 * NIST Handbook of Mathematical Functions (http://dlmf.nist.gov/18)
