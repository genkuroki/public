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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using Primes: Primes
using SymPy
@vars x

# %%
f = prod(x - k for k in 1:5) + 1

# %%
vcat([0 f.expand()], ([p sympy.factor(f, modulus=p)] for p in Primes.primes(100))...)

# %%
f = prod(x - k for k in 1:5) + 1

# %%
g = f(x+3).expand()

# %%
@vars a b

# %%
x2 = a*x+b

# %%
x3 = (x*x2).expand()(x^2=>x2).expand()

# %%
x4 = (x2^2).expand()(x^2=>x2).expand()

# %%
x5 = (x*x4).expand()(x^2=>x2).expand()

# %%
G = sympy.Poly(x5 - x + 1, x)

# %%
sympy.Poly(sympy.rem(x^5-x+1, x^2-a*x-b), x)

# %%
[(mod(a^4+3a^2*b+b^2-1, 5), mod(a^3*b+2a*b^2+1, 5)) for (a, b) in Iterators.product(1:4, 1:4)]

# %%
vcat([0 f.expand()], ([p sympy.factor(f+1, modulus=p)] for p in Primes.primes(100))...)

# %%
