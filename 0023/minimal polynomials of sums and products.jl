# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using SymPy

# %%
function minpoly(x, y, f, g, h, N)
    m = sympy.Poly(f, x).degree()
    n = sympy.Poly(g, y).degree()
    K = sympy.Poly(h, x).degree()
    L = sympy.Poly(h, y).degree()
    λ = symbols("λ[0:$N]")
    z = symbols("z")
    P = sum(λ[begin+k]*z^k for k in 0:N-1) + z^N
    p = P(z => h).expand()
    @time for k in K*N:-1:m
        p = p(x^k => x^k - x^(k-m)*f).expand()
    end
    @time for l in L*N:-1:n
        p = p(y^l => y^l - y^(l-n)*g).expand()
    end
    Q = sympy.Poly(p, x, y)
    C = Q.coeffs()
    @time sol = solve(C, λ)
    P = sum(get(sol, λ[begin+k], λ[begin+k])*z^k for k in 0:N-1) + z^N
    [
        sympy.Poly(f, x)
        sympy.Poly(g, y)
        sympy.Poly(h, x, y)
        sympy.Poly(P, z)
    ]
end

# %%
@vars a p x y
f = x^2 - a
g = y^2 - p
h = x + y
minpoly(x, y, f, g, h, 4)

# %%
@vars a p x y
f = x^2 - a
g = y^2 - p
h = x * y
minpoly(x, y, f, g, h, 2)

# %%
@vars a p x y
f = x^2 - a
g = y^3 - p
h = x + y
minpoly(x, y, f, g, h, 6)

# %%
@vars a p x y
f = x^2 - a
g = y^3 - p
h = x * y
minpoly(x, y, f, g, h, 6)

# %%
@vars a p q x y
f = x^2 - a
g = y^3 + p*y + q
h = x + y
minpoly(x, y, f, g, h, 6)

# %%
@vars a p q x y
f = x^2 - a
g = y^3 + p*y + q
h = x * y
minpoly(x, y, f, g, h, 6)

# %%
@vars a p q r x y
f = x^2 - a
g = y^4 + p*y^2 + q*y + r
h = x + y
minpoly(x, y, f, g, h, 8)

# %%
@vars a p q r x y
f = x^2 - a
g = y^4 + p*y^2 + q*y + r
h = x * y
minpoly(x, y, f, g, h, 8)

# %%
@vars a p q r s x y
f = x^2 - a
g = y^5 + p*y^3 + q*y^2 + r*y + s
h = x + y
minpoly(x, y, f, g, h, 10)

# %%
@vars a p q r s x y
f = x^2 - a
g = y^5 + p*y^3 + q*y^2 + r*y + s
h = (x * y)^2
minpoly(x, y, f, g, h, 5)

# %%
@vars a p x y
f = x^3 - a
g = y^3 - p
h = (x + y)^3
minpoly(x, y, f, g, h, 3)

# %%
@vars a p x y
f = x^3 - a
g = y^3 - p
h = x * y
minpoly(x, y, f, g, h, 3)

# %%
@vars a b p q x y
f = x^3 + a*x + b
g = y^3 + p*y + q
h = x + y
minpoly(x, y, f, g, h, 9)

# %%
@vars a b p q x y
f = x^3 + a*x + b
g = y^3 + p*y + q
h = x * y
minpoly(x, y, f, g, h, 9)

# %%
@vars a p q x y
f = x^3 - a
g = y^4 + p*y + q
h = x + y
minpoly(x, y, f, g, h, 12)

# %%
@vars a p q x y
f = x^3 - a
g = y^4 + p*y + q
h = x * y
minpoly(x, y, f, g, h, 12)

# %%
