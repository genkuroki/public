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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://twitter.com/genkuroki/status/1459220790444920833

# %%
using Distributions, StatsPlots

# %%
mixnormal = MixtureModel([Normal(0, 1), Normal(1, 2)], [1/3, 2/3])

# %%
X = rand(mixnormal, 10^5)
f(x) = (1/3)*pdf(Normal(0, 1), x) + (2/3)*pdf(Normal(1, 2), x)

histogram(X; norm=true, alpha=0.3, label="sample of mixnormal")
plot!(f; label="(1/3)N(x|0,1)+(2/3)N(x|1,2))", lw=2)

# %%
