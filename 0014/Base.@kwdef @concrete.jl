# -*- coding: utf-8 -*-
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

# %% [markdown]
# See https://github.com/jonniedie/ConcreteStructs.jl/issues/4

# %%
using ConcreteStructs

Base.@kwdef @concrete struct Problem g; y0; v0; tspan end
Problem(;
    g = 9.80665,
    y0 = 0.0,
    v0 = 30.0,
    tspan = (0.0, 8.0)
) = Problem(g, y0, v0, tspan)

prob1 = Problem()

# %%
prob2 = Problem(v0 = 100.0)

# %%
using MonteCarloMeasurements
prob3 = Problem(y0 = 0.0 ± 0.0, v0 = 30.0 ± 1.0)

# %%
