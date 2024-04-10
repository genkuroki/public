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

x = -1:0.1:1
y = @. sinpi(x) + 0.3randn()

A = x .^ (0:3)' # design matrix
betahat = A \ y # least squares estimate

scatter(x, y; label="data")
plot!(sinpi; label="true", ls=:dot)
plot!(x -> evalpoly(x, betahat); label="degree-3 polynomial regression")

# %%
