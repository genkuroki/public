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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
using SymPy
@syms x r s k

@show z = Sym(1//16)
@show y = (k^2-1)/48
expr = -(x-4y-z) - (22//3)*(x-y-z)*(x-3y-z-1) - (49//36)*(x-2y-z)*(x-2y-z-2) - (50//3)*(x-y-z)*(x-y-z-1)*(x-2y-z-2) - 4*(x-y-z)*(x-y-z-1)*(x-y-z-2)*(x-y-z-3)
@show expr = factor(expr)
@show sol = solve(expr, x)
[factor(48f + 1) for f in sol]

# %%
using SymPy
@syms h

expr = 5h + (22//9)*(-h)*(-3h-2) - (31//36)*(-2h)*(-2h-3) - (16//27)*(-(-h)*(-h-2)*(-h-4))
@show expr = factor(expr)
sol = solve(expr, h)

# %%
