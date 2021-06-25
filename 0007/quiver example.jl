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
using Plots
using Random
Random.seed!(4649373)

n = 10
x, y = randn(n+1), rand(n+1)
u, v = diff(x), diff(y)
quiver(x[1:n], y[1:n]; quiver=(u, v))

# %%
plot(x, y; arrow=true)

# %%
xs = ys = -1:0.2:1
X, Y = [x for y in ys, x in xs], [y for y in ys, x in xs]
U, V = -0.3Y, 0.3X
quiver(vec(X-U/2), vec(Y-V/2); quiver=(vec(U), vec(V)), size=(420, 400))

# %%
