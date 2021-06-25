# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
g = (k^2 for k in -3:3)

# %%
collect(g)

# %%
v = [k^2 for k in -3:3]

# %%
Tuple(g)

# %%
t = Tuple(k^2 for k in -3:3)

# %%
Set(g)

# %%
s = Set(k^2 for k in -3:3)

# %%
d = Dict(k => k^2 for k in -3:3)

# %%
:([k^2 for k in -3:3]) |> Meta.show_sexpr

# %%
:(Tuple(k^2 for k in -3:3)) |> Meta.show_sexpr

# %%
:(Set(k^2 for k in -3:3)) |> Meta.show_sexpr

# %%
:(Dict(k => k^2 for k in -3:3)) |> Meta.show_sexpr

# %%
f(x, y) = (d = x^2+y^2; d == 0 ? d : x^2*y/d)
X = Y = range(-1, 1; length=201)

z = [f(x, y) for y in Y, x in X]
Z = f.(X', Y)
z == Z

# %%
using Plots
surface(X, Y, Z; camera=(40, 70))

# %%
surface(X, Y, f; camera=(40, 70))

# %%
