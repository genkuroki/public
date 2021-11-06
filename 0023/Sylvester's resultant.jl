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

R, (a, b, c, l, m, n) = ZZ["a", "b", "c", "l", "m", "n"]
Rx, x = R["x"]

# %%
f = a*x^2 + b*x + c

# %%
g = l*x^2 + m*x + n

# %%
resultant(f, g)

# %%
