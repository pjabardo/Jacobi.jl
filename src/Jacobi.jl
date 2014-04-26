module Jacobi

# package code goes here
export jacobi, djacobi, jacobi_zeros
export zgj, zglj, zgrjm, zgrjp, wgj, wglj, wgrjm, wgrjp
export dgj, dglj, dgrjm, dgrjp
export lagrange, interp_mat



include("jac_poly.jl")
include("gauss_quad.jl")


end # module
