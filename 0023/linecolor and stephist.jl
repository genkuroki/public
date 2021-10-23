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

# %%
using Plots
gr()
X = 1100 .+ 200randn(10^6)
histogram(X; norm=true, alpha=0.3, linealpha=0.3, linecolor=:match)

# %%
stephist(X; norm=true)

# %%
