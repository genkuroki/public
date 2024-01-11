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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
x = 1:5
collect(x)

# %%
k = (0:2)'
collect(k)

# %%
A = (1:5) .^ (0:2)'

# %%
y = @. 1 + 2x + 3x^2+ 0.3rand()

# %%
betahat = A \ y

# %%
using StatsPlots

x = range(0, 2.5, 26)
f(x) = x^2 - 2x + 1
y = @. f(x) + 0.3randn()

A = x .^ (0:2)'
betahat = A \ y 

scatter(x, y; label="data", c=1)
plot!(f; label="true curve", c=1, ls=:dash)
plot!(x -> evalpoly(x, betahat); label="regression curve", c=2, lw=2)

# %%
X = 10^8:-1:1
f(x) = sin(x)/x

1 + 2sum(f.(X))
@time 1 + 2sum(f.(X))
@time 1 + 2sum(f.(X))
@time 1 + 2sum(f.(X))

# %%
1 + 2sum(f, X)
@time 1 + 2sum(f, X)
@time 1 + 2sum(f, X)
@time 1 + 2sum(f, X)

# %%
using LoopVectorization

1 + 2vsum(f, X)
@time 1 + 2vsum(f, X)
@time 1 + 2vsum(f, X)
@time 1 + 2vsum(f, X)

# %%
maximum(sin, X)

# %%
minimum(sin, X)

# %%
using Statistics

mu = mean(sin, X)

# %%
