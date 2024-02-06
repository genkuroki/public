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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

Random.seed!(1234)

x = rand(Uniform(-1, 1), 50)
u = 0.5quantile.(Normal(), abs.(x))
y = @. 2 + x + u

A = x .^ (0:1)'
@show betahat = A \ y

scatter(x, y; msc=:auto, label="data")
plot!(x -> evalpoly(x, betahat); label="regression line")
plot!(xguide="x", yguide="y")

# %%
residual = y - evalpoly.(x, (betahat,))
scatter(x, residual; msc=:auto, label="", title="residual error")
hline!([0.0]; label="")
plot!(xguide="x")

# %%
residual = y - evalpoly.(x, (betahat,))
qqnorm(residual; c=[2 1], msc=:auto, title="QQ plot of residual error")
plot!(size=(400, 400))

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

Random.seed!(1234)

t = rand(Uniform(-1, 1), 50)
u = 0.5quantile.(Normal(), abs.(t))
umin = minimum(u)
x = @. sign(t) * √(u - umin)
y = @. 2 + x + u

A = x .^ (0:1)'
@show betahat = A \ y

scatter(x, y; msc=:auto, label="data")
plot!(x -> evalpoly(x, betahat); label="regression line")
plot!(xguide="x", yguide="y")

# %%
residual = y - evalpoly.(x, (betahat,))
scatter(x, residual; msc=:auto, label="", title="residual error")
hline!([0.0]; label="")
plot!(xguide="x")

# %%
residual = y - evalpoly.(x, (betahat,))
qqnorm(residual; c=[2 1], msc=:auto, title="QQ plot of residual error")
plot!(size=(400, 400))

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

Random.seed!(1234)

t = rand(Uniform(-1, 1), 200)
u = 0.5quantile.(Normal(), abs.(t))
umin = minimum(u)
x = @. sign(t) * √(u - umin)
y = @. 2 + x + u

A = x .^ (0:1)'
@show betahat = A \ y

scatter(x, y; msc=:auto, label="data", ms=2)
plot!(x -> evalpoly(x, betahat); label="regression line")
plot!(xguide="x", yguide="y")

# %%
residual = y - evalpoly.(x, (betahat,))
scatter(x, residual; msc=:auto, label="", title="residual error", ms=2)
hline!([0.0]; label="")
plot!(xguide="x")

# %%
residual = y - evalpoly.(x, (betahat,))
qqnorm(residual; c=[2 1], msc=:auto, title="QQ plot of residual error", ms=2)
plot!(size=(400, 400))

# %%
histogram(residual; norm=true, alpha=0.3, label="", title="residual error")
plot!(fit(Normal, residual); label="")

# %%
