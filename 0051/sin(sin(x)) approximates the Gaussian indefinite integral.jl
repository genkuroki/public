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
#     display_name: Julia 1.10.7
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using QuadGK
using Plots
default(fmt=:png)

f(x) = quadgk(t -> exp(-t^2), 0, x)[1]
g(x) = sin(sin(x))

plot(f, -π/2, π/2; label="\$\\int_0^x \\exp(-t^2)\\,dt\$")
plot!(g; label="\$\\sin(\\sin(x))\$", ls=:dash)
plot!(xguide="\$x\$", guidefontsize=16)
plot!(legendfontsize=16)

# %%
