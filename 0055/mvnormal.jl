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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
using Distributions
using HCubature
using StatsPlots
default(fmt=:png)

# %%
mvn = MvNormal(zeros(2), [2 1; 1 2])
f(x) = pdf(mvn, x)
n = 10^4
X = rand(mvn, n)
scatter(X[1,:], X[2,:]; label="", msc=:auto, ms=1, ma=0.5, size=(400, 400))

# %%
c = 3.129
@show C = fill(c, 2)
@show hcubature(f, -C, C)
P = scatter(X[1,:], X[2,:]; label="", msc=:auto, ms=1, ma=0.5, size=(400, 400))
plot!(xtick=-6:6, ytick=-6:6)
plot!([-c, c, c, -c, -c], [-c, -c, c, c, -c]; label="")

# %%
d = cquantile(Normal(0, sqrt(2)), 0.05/2)
@show D = [d, 100.0]
@show hcubature(f, -D, D)
Q = scatter(X[1,:], X[2,:]; label="", msc=:auto, ms=1, ma=0.5, size=(400, 400))
plot!(xtick=-6:6, ytick=-6:6)
vline!([-d, d]; label="")

# %%
plot(P, Q; size=(400, 800), layout=(2, 1))

# %%
