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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
平均(p, X) = (sum(x -> x^p, X)/length(X))^(1/p)
算術平均(X) = sum(X)/length(X)
幾何平均(X) = prod(X)^(1/length(X))
調和平均(X) = length(X)/sum(x -> 1/x, X)

# %%
X = rand(20)
@show X
println()
@show 平均(1, X) 算術平均(X)
println()
@show 平均(1e-8, X) 幾何平均(X) 平均(-1e-8, X)
println()
@show 平均(-1, X) 調和平均(X);

# %%
