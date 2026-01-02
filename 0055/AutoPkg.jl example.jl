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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# Google Colab: https://colab.research.google.com/github/genkuroki/public/blob/main/0055/AutoPkg.jl%20example.ipynb

# %%
import Pkg
_projdeps = keys(Pkg.project().dependencies)
"AutoPkg" ∈ _projdeps || Pkg.add("AutoPkg"); using AutoPkg
AutoPkg._running_in_colab() = true # does no use Pkg.activate(; temp = true)
push!(AutoPkg._autopkg_added, _projdeps..., readdir(Sys.STDLIB)...)
@autopkg begin
    using Distributions
    using StatsPlots
end
default(fmt=:png)

# %%
dist = Gamma(4, 1)
X = rand(dist, 2^12)
histogram(X; norm=true, alpha=0.3, label="sample")
density!(X; label="kernel density estimation")
plot!(1dist; label="$dist", ls=:dash)

# %%
