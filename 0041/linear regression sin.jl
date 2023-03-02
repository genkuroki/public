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
#     display_name: Julia 1.9.0-beta4
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
n = 10^6
x = rand(Uniform(0, 10), n)
y = @. sinpi(x) + rand(Normal())
stephist(y; norm=true, label="")
plot!(fit(Normal, y); label="")

# %%
n = 10^3
x = rand(Uniform(0, 4), n)
y = @. 3 + 5x + sinpi(x) + rand(Normal())

scatter(x, y; msc=:auto, alpha=0.5, ms=3, label="")
#plot!(x -> 1 + x + sinpi(x); label="", lw=1.5, c=:blue)

# %%
X = [ones(n) x]
@show betahat0, betahat1 = X \ y

scatter(x, y; msc=:auto, alpha=0.5, ms=3, label="")
plot!(x -> betahat0 + betahat1*x; label="", lw=1.5)
#plot!(x -> 1 + x + sinpi(x); label="", lw=1.5, c=:blue)

# %%
yhat = @. betahat0 + betahat1*x
uhat = y - yhat
normal = fit(Normal, uhat)

stephist(uhat; norm=true, label="")
plot!(normal; label="")

# %%
scatter(x, uhat; msc=:auto, alpha=0.5, ms=3, label="")
hline!([0]; label="", lw=1.5)
#plot!(sinpi; label="", lw=1.5, c=:blue)

# %%
