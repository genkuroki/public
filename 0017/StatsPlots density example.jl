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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
using Plots
using Distributions

# %%
d = MixtureModel([Poisson(100), Poisson(250)], [0.95, 0.05])
n = 10^4
sample = sum.([rand(d, 20) for _ in 1:n])
histogram(sample)

# %%
using StatsPlots
x=rand(1:6, 10^5)
density(x; label="default")
density!(x; label="bandwidth=1", bandwidth=1)

# %%
using StatsPlots
x = rand(1:6, 10^5)
histogram(x; norm=true, alpha=0.3, bin=0.5:6.5, label="histogram of true sample")
density!(x; label="kde with bandwidth=1", bandwidth=1, lw=2)
plot!(; ylim=(-0.005, 0.22), xtick=-10:10)

# %%
