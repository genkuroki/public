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
@vars x y z

# %% [markdown]
# ## $\theta = \sqrt{2} + \sqrt[3]{3}$ の $\mathbb{Q}$ 上での最小多項式

# %%
a = √Sym(2)

# %%
b = Sym(3)^(1/Sym(3))

# %%
θ = a + b

# %%
sympy.minimal_polynomial(θ, x)

# %%
ω = (-1 + √Sym(-3))/2

# %%
ω^2 + ω + 1 |> expand

# %%
prod((x - (a + ω^k * b))*(x - (-a + ω^k * b)) for k in 0:2).expand()

# %%
f = sympy.minimal_polynomial(a, x)

# %%
g = sympy.minimal_polynomial(b, x)

# %%
h = f(z - x).expand()

# %%
r = g - (x + 2z)*h |> expand

# %%
Denom, Numer = sympy.Poly(r, x).coeffs() |> C -> [C[1], -C[2]]

# %%
Numer / Denom

# %%
denom = Denom(z => θ).expand().simplify()

# %%
numer = Numer(z => θ).expand().simplify()

# %%
c = numer / denom |> simplify

# %%
c / b |> simplify

# %% [markdown]
# ## $\omega = (-1 + \sqrt{-3})/2$ のときの $\theta = \omega + \sqrt[3]{2}$ の $\mathbb{Q}$ 上での最小多項式

# %%
ω = (-1 + √Sym(-3))/2

# %%
ω^2 + ω + 1 |> expand

# %%
b = Sym(2)^(1/Sym(3))

# %%
θ = ω + b

# %%
sympy.minimal_polynomial(θ, x)

# %%
prod((x - (ω + ω^k * a))*(x - (ω^2 + ω^k * a)) for k in 0:2).expand()

# %%
f = sympy.minimal_polynomial(ω, x)

# %%
g = sympy.minimal_polynomial(b, x)

# %%
h = f(z - x).expand()

# %%
r = g - (x + 2z + 1)*h |> expand

# %%
Denom, Numer = sympy.Poly(r, x).coeffs() |> C -> [C[1], -C[2]]

# %%
Numer / Denom

# %%
denom = Denom(z => θ) |> expand |> simplify

# %%
numer = Numer(z => θ) |> expand |> simplify

# %%
c = numer/denom |> simplify

# %%
c / b |> simplify

# %% [markdown]
# ## $\theta = \sqrt{-3} + \sqrt[3]{2}$ の $\mathbb{Q}$ 上での最小多項式

# %%
a = √Sym(-3)

# %%
b = Sym(2)^(1/Sym(3))

# %%
θ = a + b

# %%
sympy.minimal_polynomial(θ, x)

# %%
ω = (-1 + √Sym(-3))/2

# %%
ω^2 + ω + 1 |> expand

# %%
prod((x - (a + ω^k * b))*(x - (-a + ω^k * b)) for k in 0:2).expand()

# %%
f = sympy.minimal_polynomial(a, x)

# %%
g = sympy.minimal_polynomial(b, x)

# %%
h = f(z - x).expand()

# %%
r = g - (x + 2z)*h |> expand

# %%
Denom, Numer = sympy.Poly(r, x).coeffs() |> C -> [C[1], -C[2]]

# %%
Numer / Denom

# %%
denom = Denom(z => θ) |> expand |> simplify

# %%
numer = Numer(z => θ) |> expand |> simplify

# %%
c = numer/denom |> simplify

# %%
c / b |> simplify

# %% [markdown]
# ## $\theta = \sqrt{2} + \sqrt[5]{3}$ の $\mathbb{Q}$ 上での最小多項式

# %%
a = √Sym(2)

# %%
b = Sym(3)^(1/Sym(5))

# %%
θ = a + b

# %%
sympy.minimal_polynomial(θ, x)

# %%
ω = cos(2PI/5) + im*sin(2PI/5)

# %%
ω^4 + ω^3 + ω^2 + ω + 1 |> expand

# %%
prod((x - (a + ω^k * b))*(x - (-a + ω^k * b)) for k in 0:4).expand()

# %%
f = sympy.minimal_polynomial(a, x)

# %%
g = sympy.minimal_polynomial(b, x)

# %%
h = f(z - x).expand()

# %%
r = g - (x^3 + 2z*x^2 + (3z^2+2)*x + 4z^3+8z)*h |> expand

# %%
Denom, Numer = sympy.Poly(r, x).coeffs() |> C -> [C[1], -C[2]]

# %%
Numer / Denom

# %%
denom = Denom(z => θ) |> expand |> simplify

# %%
numer = Numer(z => θ) |> expand |> simplify

# %%
c = numer/denom |> simplify

# %%
c / b |> simplify

# %%
