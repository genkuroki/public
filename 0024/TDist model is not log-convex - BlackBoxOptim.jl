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
#     display_name: Julia 1.6.4
#     language: julia
#     name: julia-1.6
# ---

# %%
using Distributions
using Plots
using BlackBoxOptim

# model
Distributions.TDist(μ, ρ, ν) = LocationScale(μ, ρ, TDist(ν))

# test data
X = [-0.01, 0.01, 1.0];

# %%
using Random
Random.seed!(123456789)

# %%
SearchRange = [(-0.1, 0.4), (-15.0, 10.0), (-2.5, 8.0)]
TraceMode = :compact
TraceInterval = 0.001
MaxSteps = 50000
@show SearchRange
o = bboptimize(x -> -loglikelihood(TDist(x[1], 10^x[2], 10^x[3]), X); SearchRange, TraceMode, TraceInterval, MaxSteps);

# %%
SearchRange = [(-0.5, 1.0), (-15.0, 10.0), (-2.5, 8.0)]
TraceInterval = 0.001
MaxSteps = 50000
@show SearchRange
o = bboptimize(x -> -loglikelihood(TDist(x[1], 10^x[2], 10^x[3]), X); SearchRange, TraceMode, TraceInterval, MaxSteps);

# %%
SearchRange = [(-3.0, 4.0), (-15.0, 10.0), (-2.5, 8.0)]
TraceInterval = 0.001
MaxSteps = 50000
@show SearchRange
o = bboptimize(x -> -loglikelihood(TDist(x[1], 10^x[2], 10^x[3]), X); SearchRange, TraceMode, TraceInterval, MaxSteps);

# %%
