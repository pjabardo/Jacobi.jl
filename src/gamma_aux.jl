function rgammanm(z, n, m=0)
    if n==m
        return one(z)
    elseif n < m
        return one(z) / rgammanm(z, m, n)
    end
        
    rg = one(z)
    for i = (n-1):-1:m
        rg = rg * (z + i)
    end

    return rg
end

function rgammaxy(x, y, n)
    rg = x/y * gamma(x)/gamma(y)

    for i = 1:(n-1)
        rg = rg * (x+i)/(y+i)
    end
    return rg
end


function rgammaxnym(x, n, y, m)

    if x == y
        return rgammanm(x, n, m)
    elseif n == m
        return rgammaxy(x, y, n)
    elseif n < m
        return one(x)/rgammaxnym(y, m, x, n)
    end

    rg = x/y * gamma(x)/gamma(y)

    for i = (m-1):-1:1
        rg = rg * (x+i)/(y+i)
    end

    for i = (n-1):-1:m
        rg = rg * (x+i)
    end

    return rg
end

