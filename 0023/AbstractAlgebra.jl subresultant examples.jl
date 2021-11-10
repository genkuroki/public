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

using AbstractAlgebra
AbstractAlgebra.PrettyPrinting.set_html_as_latex(true)

eqstr(x) = "\\ \\text{is " * string(x) * "}"
eqstr(x::SetElem) = "=" * sprint(x) do io, x; show(io, "text/latex", x) end

function dispeq(lhs, rhs)
    latexstr = "\$\\displaystyle $lhs " * eqstr(rhs) * "\$"
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
    C = coefficients(f)
    sum(a*z^(i-1)*x^(m-i+1) for (i, a) in enumerate(C))
end

function val(h::FracElem, θ, L)
    n = L(numerator(h)(θ))
    d = L(denominator(h)(θ))
    n//d
end

function testval(h::FracElem, θ, β, L)
    n = L(numerator(h)(θ))
    d = L(denominator(h)(θ))
    n == β*d
end

function calcallresults(z, x, f, g)
    f_plus = (-1)^degree(f)*f(z - x)
    F_plus = resultant(f_plus, g)
    β_plus = commonroot(g, f_plus)
    
    f_mult = (-1)^degree(f)*revpoly(z, x, f)
    F_mult = resultant(f_mult, g)
    β_mult = commonroot(g, f_mult)
    
    K = FractionField(base_ring(z))
    B1, α = K["α"]
    L1 = ResidueField(B1, numerator(f(z))(α))
    B2, β = L1["β"]
    L = ResidueField(B2, numerator(g(z))(β))
    
    val_F_plus = val(F_plus, α + β, L)
    test_β_plus = testval(β_plus, α + β, β, L)
    val_F_mult = val(F_mult, α * β, L)
    test_β_mult = testval(β_mult, α * β, β, L)
    
    S_plus = subresultant((-1)^degree(f)*f(z-x), g, 1)
    rootS_plus = rootdeg1(S_plus)
    
    F_plus, β_plus, F_mult, β_mult, val_F_plus, test_β_plus, val_F_mult, test_β_mult, S_plus, rootS_plus
end

function dispallresults(z, x, f, g)
    F_plus, β_plus, F_mult, β_mult, val_F_plus, test_β_plus, val_F_mult, test_β_mult, S_plus, rootS_plus =
        @time calcallresults(z, x, f, g)
    flush(stdout)
    dispeq("F_\\alpha(x)", f)
    dispeq("F_\\beta(x)", g)
    dispeq("R_{\\alpha + \\beta}(z)", F_plus)
    dispeq("R_{\\alpha + \\beta}(\\alpha + \\beta)", val_F_plus)
    dispeq("\\beta_\\mathrm{plus}(z)", β_plus)
    dispeq("\\beta_\\mathrm{plus}(\\alpha + \\beta) = \\beta", test_β_plus)
    dispeq("R_{\\alpha\\beta}(z)", F_mult)
    dispeq("R_{\\alpha\\beta}(\\alpha\\beta)", val_F_mult)
    dispeq("\\beta_\\mathrm{mult}(z)", β_mult)
    dispeq("\\beta_\\mathrm{mult}(\\alpha\\beta) = \\beta", test_β_mult)
    dispeq("\\text{1-subresultant}", S_plus)
    dispeq("\\text{root of 1-subresultant}", rootS_plus)
end

safecoeff(f, k) = 0 ≤ k ≤ degree(f) ? coeff(f, k) : zero(base_ring(f))

function subresultant_matrix(p::PolyElem{T}, q::PolyElem{T}, k) where T <: RingElement
    check_parent(p, q)
    R = parent(p)
    m = degree(p)
    n = degree(q)
    if length(p) == 0 || length(q) == k || m + n < 2k
       return zero_matrix(R, 0, 0)
    end
    M = zero_matrix(R, m + n - 2k, m + n - 2k)
    x = gen(R)
    for i in 1:n-k
        for j in 1:m+n-2k-1
           M[i, j] = safecoeff(p, m + (i-1) - (j-1))
        end
        M[i, end] = x^(n-k-i)*p
    end
    for i in 1:m-k
        for j in 1:m+n-2k-1
           M[n-k+i, j] = safecoeff(q, n + (i-1) - (j-1))
        end
        M[n-k+i, end] = x^(m-k-i)*q
    end
    return M
end

subresultant(p, q, k) = det(subresultant_matrix(p, q, k))

# %%
R, (a, b, c, p, q, r, s) = ZZ["a", "b", "c", "p", "q", "r", "s"]
K = FractionField(R)
Rz, z = R["z"]
Kz = FractionField(Rz)
Rx, x = Kz["x"]

# %%
f = x^3 + a*x^2 + b*x + c
g = x^4 + p*x^3 + q*x^2 + r*x + s
subresultant_matrix(f, g, 0)

# %%
S0 = sylvester_matrix(f, g)

# %%
subresultant(f, g, 0) == resultant(f, g) == det(S0)

# %% tags=[]
subresultant_matrix(f, g, 1)

# %%
S1 = Rx[
    1 a b c 0
    0 1 a b c*x
    0 0 1 a b*x+c
    1 p q r s*x
    0 1 p q r*x+s
]

# %%
subresultant(f, g, 1)

# %% tags=[]
subresultant(f, g, 1) == det(S1)

# %%
S2 = Rx[
    1 a b*x^2+c*x
    0 1 a*x^2+b*x+c
    1 p q*x^2+r*x+s
]

# %%
subresultant(f, g, 2)

# %%
subresultant(f, g, 2) == det(S2)

# %%
f = x^3 - a
g = x^4 - p
dispallresults(z, x, f, g)

# %%
subresultant_matrix(-f(z-x), g, 1)

# %%
s1 = Rx[
    1 -3z 3z^2 -z^3+a 0
    0 1   -3z  3z^2   (-z^3+a)*x
    0 0   1    -3z    3z^2*x-z^3+a
    1 0   0    0      -p*x
    0 1   0    0      -p
]

# %%
S_plus = subresultant(-f(z-x), g, 1)

# %%
subresultant(-f(z-x), g, 1) == det(s1)

# %%
rootS_plus = rootdeg1(S_plus)
dispeq("\\text{root of 1-subresultant}", rootS_plus)

# %% tags=[]
f = x^3 - a
g = x^3 + p*x + q
dispallresults(z, x, f, g)

# %% tags=[]
f = x^3 - a
g = x^3 + p*x + q
dispallresults(z, x, g, f)

# %%
h1 = g(z-x) + f

# %%
(3z)^2*f - 3z*x*h1

# %%
h2 = (3z)^2*f - (3z*x + 3z^2 + p)*h1

# %%
subresultant(f, g(z-x), 1)

# %%
f = x^2 - a
g = x^2 - p
dispallresults(z, x, f, g)

# %%
f = x^2 - a
g = x^3 - p
dispallresults(z, x, f, g)

# %%
f = x^2 - a
g = x^4 - p
dispallresults(z, x, f, g)

# %%
f = x^2 - a
g = x^5 - p
dispallresults(z, x, f, g)

# %%
f = x^3 - a
g = x^3 - p
dispallresults(z, x, f, g)

# %%
f = x^3 - a
g = x^4 - p
dispallresults(z, x, f, g)

# %%
f = x^3 - a
g = x^5 - p
dispallresults(z, x, f, g)

# %%
f = x^3 + a*x + b
g = x^3 + p*x + q
dispallresults(z, x, f, g)

# %%
f = x^2 - a
g = x^3 + p*x^2 + q*x + r
dispallresults(z, x, f, g)

# %%
f = x^3 + a*x + b
g = x^4 + p*x^2 + q*x + r
dispallresults(z, x, f, g)

# %%
