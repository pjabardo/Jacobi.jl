module Jacobi

# package code goes here
export jacobi, djacobi, jacobi_zeros
export zgj, zglj, zgrjm, zgrjp, wgj, wglj, wgrjm, wgrjp
export dgj, dglj, dgrjm, dgrjp, lgj, lglj, lgrjm, lgrjp
export igj, iglj, igrjm, igrjp



include("jac_poly.jl")
include("gauss_quad.jl")


end # module
