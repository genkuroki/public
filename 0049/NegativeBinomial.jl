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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

r, p = 3, 0.2
@show r p
negbin = NegativeBinomial(r, p) + r/2
@show μ = mean(negbin);
@show σ² = var(negbin);
@show gam = Gamma(μ^2/σ², σ²/μ);
bar(negbin; alpha=0.3, label="NegativeBinomial(r, p) + r/2")
plot!(gam; label="Gamma(μ^2/σ², σ²/μ)")
title!("r=$r, p=$p")

# %%
?NegativeBinomial

# %%
