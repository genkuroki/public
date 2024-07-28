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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# * https://x.com/penta_math/status/1817130895893287136
# * https://x.com/dannchu/status/1817456970514993179

# %%
using Distributions
using StatsPlots
default(fmt=:png)

X = rand(10^6)
Y = @. -log(X)
stephist(Y; norm=true, label="Y=-log(X) where X~Uniform(0,1)")
plot!(Exponential(); label="Exponential(1)", ls=:dash, lw=2)
plot!(size=(400, 250), xlim=(-0.2, 6))

# %%
function sim(; dist=Uniform(0, 10), f=prod, n=5, a=1000, L=10^8)
    count(f(rand(dist) for _ in 1:n) ≤ a for _ in 1:L) / L
end

sim(dist=Uniform(0, 10), f=prod, n=5, a=1000)

# %%
sim(dist=Uniform(0, 1), f=prod, n=5, a=0.01)

# %%
1 - sim(dist=Exponential(), f=sum, n=5, a=-log(0.01))

# %%
1 - sim(dist=Gamma(5, 1), f=sum, n=1, a=-log(0.01))

# %%
ccdf(Gamma(5, 1), -log(0.01))

# %%
