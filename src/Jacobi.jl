module Jacobi

# package code goes here
export jacobi, djacobi, jacobi_zeros, jacobi_zeros!
export zgj, zglj, zgrjm, zgrjp, wgj, wglj, wgrjm, wgrjp
export dgj, dglj, dgrjm, dgrjp
export lagrange, interp_mat
export chebyshev, chebyshev2, dchebyshev, dchebyshev2
export chebyshev_zeros!, chebyshev_zeros, chebyshev2_zeros!, chebyshev2_zeros
export legendre, dlegendre, legendre_zeros
export chebyshev!, chebyshev2!, dchebyshev!, dchebyshev2!
export legendre!, jacobi!, djacobi!

include("jac_poly.jl")
include("chebyshev.jl")
include("legendre.jl")
include("gauss_quad.jl")


end # module
