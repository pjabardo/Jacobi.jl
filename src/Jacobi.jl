module Jacobi

# package code goes here
export jacobi, djacobi, jacobi_zeros, jacobi_zeros!
export zgj, zglj, zgrjm, zgrjp, wgj, wglj, wgrjm, wgrjp
export dgj, dglj, dgrjm, dgrjp
export lagrange, lagrange!, interp_mat
export chebyshev, chebyshev2, dchebyshev, dchebyshev2
export chebyshev_zeros!, chebyshev_zeros, chebyshev2_zeros!, chebyshev2_zeros
export legendre, dlegendre, legendre_zeros
export chebyshev!, chebyshev2!, dchebyshev!, dchebyshev2!
export legendre!, jacobi!, djacobi!
export poly_chebyshev, poly_chebyshev2, poly_dchebyshev, poly_dchebyshev2
export poly_legendre, poly_dlegendre, poly_jacobi, poly_djacobi

export QUADRATURE_TYPE, GJ, GLJ, GRJM, GRJP
export qzeros, qweights, qdiff, Quadrature
export qalpha, qbeta



include("jac_poly.jl")
include("chebyshev.jl")
include("legendre.jl")
include("gauss_quad.jl")
include("polynomials.jl")

end # module
