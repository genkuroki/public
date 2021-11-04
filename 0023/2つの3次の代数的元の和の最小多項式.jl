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
Poly = sympy.Poly
@vars a b p q x α β

# %%
λ = symbols("λ[0:9]")

# %%
f = x^3 + a*x + b
Poly(f, x)

# %%
g = x^3 + p*x + q
Poly(g, x)

# %%
H = sum(λ[begin+k]*x^k for k in 0:8) + x^9
Poly(H, x)

# %%
h = H(x => α + β).expand()
for k in 6:-1:0
    h = h(α^(k+3) => α^(k+3) - α^k * f(x=>α), β^(k+3) => β^(k+3) - β^k * g(x=>β)).expand()
end
P = Poly(h, α, β)

# %%
C = P.coeffs()

# %%
length(C)

# %%
sol = @time solve(C, λ)

# %% [markdown]
# Assume that
#
# * the minimal polynomial of $\alpha$ is equal to $x^3 + ax + b$,
# * the minimal polynomial of $\beta$ is equal to $x^3 + px + q$, and
# * the minimal polynomial of $\alpha + \beta$ is of degree $9$.
#
# Then the minimal polynomial of $\alpha + \beta$ is equal to $x^{9} + \left(3 a + 3 p\right) x^{7} + \left(3 b + 3 q\right) x^{6} + \left(3 a^{2} + 3 a p + 3 p^{2}\right) x^{5} + \left(6 a b - 3 a q - 3 b p + 6 p q\right) x^{4} + \left(a^{3} + a^{2} p + a p^{2} + 3 b^{2} - 21 b q + p^{3} + 3 q^{2}\right) x^{3} + \left(3 a^{2} b + 3 a^{2} q - 6 a b p - 6 a p q + 3 b p^{2} + 3 p^{2} q\right) x^{2} + \left(a^{3} p - 2 a^{2} p^{2} + 3 a b^{2} - 3 a b q + a p^{3} - 6 a q^{2} - 6 b^{2} p - 3 b p q + 3 p q^{2}\right) x + a^{3} q + a^{2} b p - 2 a^{2} p q - 2 a b p^{2} + a p^{2} q + b^{3} + 3 b^{2} q + b p^{3} + 3 b q^{2} + q^{3}$.

# %%
h = sum(sol[λ[begin+k]]*x^k for k in 0:8) + x^9
Poly(h, x)

# %% [markdown]
# Assume that
#
# * the minimal polynomial of $\alpha$ is equal to $x^3 + ax + b$,
# * the minimal polynomial of $\beta$ is equal to $x^3 - p $, and
# * the minimal polynomial of $\alpha + \beta$ is of degree $9$.
#
# Then the minimal polynomial of $\alpha + \beta$ is equal to $x^{9} + 3 a x^{7} + \left(3 b - 3 p\right) x^{6} + 3 a^{2} x^{5} + \left(6 a b + 3 a p\right) x^{4} + \left(a^{3} + 3 b^{2} + 21 b p + 3 p^{2}\right) x^{3} + \left(3 a^{2} b - 3 a^{2} p\right) x^{2} + \left(3 a b^{2} + 3 a b p - 6 a p^{2}\right) x -  a^{3} p + b^{3} - 3 b^{2} p + 3 b p^{2} - p^{3}$

# %%
h1 = h(p=>0)(q=>-p)
Poly(h1, x)

# %% [markdown]
# Assume that
#
# * the minimal polynomial of $\alpha$ is equal to $x^3 - a$,
# * the minimal polynomial of $\beta$ is equal to $x^3 - b$, and
# * the minimal polynomial of $\alpha + \beta$ is of degree $9$.
#
# Then the minimal polynomial of $\alpha + \beta$ is equal to $x^{9} + \left(- 3 a - 3 b\right) x^{6} + \left(3 a^{2} - 21 a b + 3 b^{2}\right) x^{3} -  a^{3} - 3 a^{2} b - 3 a b^{2} - b^{3}$.

# %%
h2 = h(a=>0, p=>0)(b=>-a)(q=>-b)
Poly(h2, x)

# %% [markdown]
# The minimal polynomial of $\sqrt[3]{2} + \sqrt[3]{3}$ is equal to $x^{9} - 15 x^{6} - 87 x^{3} - 125$.

# %%
h3 = h2(a=>2, b=>3)
Poly(h3, x)

# %%
