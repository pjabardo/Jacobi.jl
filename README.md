# Jacobi

A library that implements Jacobi polynomials and Gauss 
quadrature related operations.

## Jacobi polynomials

Jacobi polynomials  are implemented in function 
`jacobi(x, n, a, b)` where `x` is where you want to calculate the 
polynomial of degree `n` with weights `a` and `b`.

To calculate the derivative of jacobi polynomials, use function 
`djacobi(x, n, a, b)`. 

When calculating the weights of Gaussian qudrature, it is necessary to determine
the zeros of Jacobi polynomials. The function `jacobi_zeros(m, a, b)` calculates the 
zeros of a Jacobi polynomial of degree `m` with weights `a` and `b`.


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




[![Build Status](https://travis-ci.org/pjabardo/Jacobi.jl.png)](https://travis-ci.org/pjabardo/Jacobi.jl)
