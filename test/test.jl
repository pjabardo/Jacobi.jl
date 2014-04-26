include("../src/Jacobi.jl")
using Jacobi

function test_integrate(zfun, wfun, ifun, T=Float64, qmax=20, qmin=3)
    a = zero(T)
    b = zero(T)
    ixe = ifun(one(T)) - ifun(-one(T))
    r = [qmin:qmax]
    n = length(r)
    err = zeros(T, n)
    for i = 1:n
        q = r[i]
        z = zfun(q, a, b)
        w = wfun(z, a, b)
        ix = sum(w .* fun(z))
        err[i] = ix - ixe
    end

    err
end


function test_deriv(zfun, Dfun, fun, dfun, a=0, b=0, T=Float64, qmax=20, qmin=3)
    
    a = convert(T, a)
    b = convert(T, b)

    r = [qmin:qmax]
    n = length(r)

    err = zeros(T, n)

    for i = 1:n
        q = r[i]
        z = zfun(q, a, b)
        D = Dfun(z, a, b)
        dxe = dfun(z)
        x = fun(z)
        dx = D * x
        err[i] = maximum(abs(dx-dxe))
    end
    err
end


    
function test_interp(zfun, Ifun, zi, fun, a=0, b=0, T=Float64, qmax=20, qmin=3)
    a = convert(T, a)
    b = convert(T, b)

    r = [qmin:qmax]
    n = length(r)

    err = zeros(T, n)
    xei = fun(zi)
    np = length(zi)

    for i = 1:n
        q = r[i]
        z = zfun(q, a, b)
        I = Ifun(zi, z, a, b)
        fi = fun(z)
        xi = I * fi
        err[i] = maximum(abs(xi - xei))
    end
    err
end


fun(x) = x.*sin(2*pi*x)
ifun(x) = ( sin(2*pi*x) - 2*pi*x .* cos(2*pi*x) ) / (4*pi*pi)
dfun(x) = sin(2*pi*x) + 2*pi*x .* cos(2*pi*x)

function runtests(a=0, b=0, T=Float64, qmax=30, qmin=3)



    xi = linspace(-1, 1, 11)
    yi = fun(xi)

    zfun = [zgj, zglj, zgrjm, zgrjp]
    wfun = [wgj, wglj, wgrjm, wgrjp]
    Dfun = [dgj, dglj, dgrjm, dgrjp]
    Ifun = [igj, iglj, igrjm, igrjp]


    Q = [qmin:qmax]
    nq = length(qmin:qmax)

    nquad = length(zfun)
    ntests = 3
    err = zeros(nq, nquad, ntests)

    for i = 1:nquad
        err[:,i,1] = test_integrate(zfun[i], wfun[i], ifun, T, qmax, qmin)
        err[:,i,2] = test_deriv(zfun[i], Dfun[i], fun, dfun, a, b, T, qmax, qmin)
        err[:,i,3] = test_interp(zfun[i], Ifun[i], xi, fun, a, b, T, qmax, qmin)
    end

    return err
end




        

    