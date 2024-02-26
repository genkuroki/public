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
#     display_name: Julia 1.10.1
#     language: julia
#     name: julia-1.10
# ---

# %%
# https://x.com/fermatslibrary/status/1761482490122748355?s=61&t=_KnHkB3gSNKRbi3Ce1GncA

using Plots
default(fmt=:png)

f(x) = 16(π - x)*x/(5π^2 - 4(π - x)*x)

P = plot(f, -0.5, π+0.5; label="f(x)")
plot!(sin; label="sin(x)", ls=:dash)
plot!(legend=:topright, ytick=-1:0.1:1)

Q = plot(x -> f(x) - sin(x), 0, π; label="f(x) - sin(x)")
plot!(legend=:top)

plot(P, Q; size=(600, 800), layout=(2, 1))

# %%
