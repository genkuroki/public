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
using StatsPlots
using Distributions

# %%
n = 10^3
bin = Binomial(n, 5/n)
k = 0:15
plot(k, pdf.(bin, k))
plot!(k, pdf.(Normal(mean(bin), std(bin)), k))
plot!(k, pdf.(Poisson(5), k))

# %%
