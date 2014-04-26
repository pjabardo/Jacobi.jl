module Jacobi

# package code goes here
export jacobi, djacobi, jacobi_zeros
export zgj, zglj, zgrjm, zgrjp, wgj, wglj, wgrjm, wgrjp
export dgj, dglj, dgrjm, dgrjp, lgj, lglj, lgrjm, lgrjp
export igj, iglj, igrjm, igrjp
export rgammanm, rgammaxy, rgammaxnym



include("jac_poly.jl")
include("gauss_quad.jl")
include("gamma_aux.jl")


end # module
