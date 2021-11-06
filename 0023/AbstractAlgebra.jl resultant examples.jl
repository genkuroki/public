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
using AbstractAlgebra
AbstractAlgebra.PrettyPrinting.set_html_as_latex(true)

function dispeq(lhs, rhs)
    str = sprint(rhs) do io, x; show(io, "text/latex", x) end
    latexstr = "\$\\displaystyle $lhs =" * str * "\$"
    display("text/html", latexstr)
end

function euclidean(f, g, stop = 1)
    while true
        r = mod(f, g)
        iszero(r) && return g
        degree(r) ≤ stop && return r
        f, g = g, r
    end
end

function rootdeg1(f)
    C = coefficients(f) |> collect
    -first(C)//last(C)
end

function commonroot(f, g)
    m = euclidean(f, g, 1)
    rootdeg1(m)
end

function revpoly(z, x, f)
    m = degree(f)
    C = coefficients(f) |> collect
    sum(a*z^(i-1)*x^(m-i+1) for (i, a) in enumerate(C))
end

function calcres(z, x, f, g)
    f_plus = (-1)^degree(f)*f(z - x)
    F_plus = resultant(f_plus, g)
    β_plus = commonroot(g, f_plus)
    
    f_mult = (-1)^degree(f)*revpoly(z, x, f)
    F_mult = resultant(f_mult, g)
    β_mult = commonroot(g, f_mult)
    
    F_plus, β_plus, F_mult, β_mult
end

function dispresult(z, x, f, g)
    F_plus, β_plus, F_mult, β_mult = @time calcres(z, x, f, g)
    flush(stdout)
    dispeq("F_\\alpha(x)", f)
    dispeq("F_\\beta(x)", g)
    dispeq("R_{\\alpha + \\beta}(z)", F_plus)
    dispeq("\\beta_\\mathrm{plus}(z)", β_plus)
    dispeq("R_{\\alpha\\beta}(z)", F_mult)
    dispeq("\\beta_\\mathrm{mult}(z)", β_mult)
end

# %%
R, (a, b, c, p, q, r, s) = ZZ["a", "b", "c", "p", "q", "r", "s"]
Rz, z = R["z"]
Kz = FractionField(Rz)
Rx, x = Kz["x"]

# %%
f = x^2 - a
g = x^2 - p
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^3 - p
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^4 - p
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^5 - p
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^6 - p
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^7 - p
dispresult(z, x, f, g)

# %%
f = x^3 - a
g = x^4 - p
dispresult(z, x, f, g)

# %%
f = x^3 - a
g = x^5 - p
dispresult(z, x, f, g)

# %%
f = x^3 - a
g = x^6 - p
dispresult(z, x, f, g)

# %%
f = x^3 - a
g = x^7 - p
dispresult(z, x, f, g)

# %%
f = x^4 - a
g = x^4 - p
dispresult(z, x, f, g)

# %%
f = x^4 - a
g = x^5 - p
dispresult(z, x, f, g)

# %%
f = x^5 - a
g = x^5 - p
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^3 + p*x + q
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^4 + p*x^2 + q*x + r
dispresult(z, x, f, g)

# %%
f = x^2 - a
g = x^5 + p*x^3 + q*x^2 + r*x + s
dispresult(z, x, f, g)

# %%
f = x^3 - a
g = x^3 + p*x + q
dispresult(z, x, f, g)

# %%
f = x^3 - a
g = x^4 + p*x^2 + q*x + r
dispresult(z, x, f, g)

# %%
f = x^3 - a
g = x^5 + p*x^3 + q*x^2 + r*x + s
dispresult(z, x, f, g)

# %%
f = x^3 + a*x + b
g = x^3 + p*x + q
dispresult(z, x, f, g)

# %%
f = x^3 + a*x + b
g = x^4 + p*x^2 + q*x + r
dispresult(z, x, f, g)

# %%
