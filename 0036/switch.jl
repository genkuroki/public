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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %%
using Plots

switch(f, args...; kwargs...) = f(args...; kwargs...)

F(x) = switch(x) do x
    x < 0 && return 0.0
    x ≤ 1 && return 1.0
    return 2.0
end

plot(F, -1, 2; label="y = F(x)", legend=:topleft)

# %%
