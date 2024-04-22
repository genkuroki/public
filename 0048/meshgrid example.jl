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
using Plots
default(fmt=:png)

f(x, y) = cos(x*y)
x = range(-4, 4, 81)
y = range(-2.5, 2.5, 51)
xx, yy = [x for y in y, x in x], [y for y in y, x in x]
zz = f.(xx, yy)
xx, yy, zz

scatter(xx, yy; marker_z=zz, c=:bwr, ms=2, msw=0, label="")

# %%
using Distributions
using Plots
default(fmt=:png)

f(x, y) = cos(x*y)
n = 2^13
xx = rand(Normal(0, 4), n)
yy = rand(Normal(0, 4), n)
zz = f.(xx, yy)
xx, yy, zz

scatter(xx, yy; marker_z=zz, c=:bwr, ms=3, msw=0, label="", alpha=0.5)
plot!(xlim=(-4, 4), ylim=(-2.5, 2.5))

# %%
using Plots
default(fmt=:png)

f(x, y) = cos(x*y)
x = range(-4, 4, 81)
y = range(-2.5, 2.5, 51)
xx, yy = [x for y in y, x in x], [y for y in y, x in x]
zz = f.(xx, yy)

X = reshape(xx, 51, 81)[1,:]
Y = reshape(yy, 51, 81)[:,1]
Z = reshape(zz, 51, 81)
heatmap(X, Y, Z; c=:bwr)

# %%
