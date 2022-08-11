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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# # HDI - highest density interval
#
# * 黒木玄
# * 2022-08-06～2022-08-11
# * License: https://opensource.org/licenses/MIT
#
# 2022-08-11: [古いバージョン](https://github.com/genkuroki/public/blob/7ebc143b3280171fa81f5cb6a28502bc80f4bd1d/0034/HDI%20-%20highest%20density%20interval.ipynb)の `hdi(dist, α)` (Roots.jlを使うバージョン)はかなり不安定.  Optim.jl を使うよりシンプルな方法に切り替えた.  

# %%
using Distributions
using Optim
using StatsPlots
default(fmt=:png, size=(400, 250), titlefontsize=10)

# %%
"""
    hdi(dist, α = 0.05; alg = Brent())

returns the 100(1 - `α`)% highest density interval (HDI) of the distribution `dist`.

Assumption: `dist` is unimodal.
"""
function hdi(dist, α = 0.05; alg = Brent())
    o = optimize(0, α, alg) do p
        a, b = quantile.(dist, (p, p + (1 - α)))
        b - a
    end
    p = o.minimizer
    q = p + (1 - α)
    quantile.(dist, (p, q))
end

@doc hdi

# %%
function plot_hdi(dist, α = 0.05; alg=Brent(), kwargs...)
    @show α
    @show a, b = hdi(dist, α; alg)
    plot(dist; label="")
    vline!([a, b]; label="Optim", ls=:dash)
    title!("HDI of $dist")
    plot!(; kwargs...)
end

# %%
plot_hdi(Normal())

# %%
plot_hdi(Gamma(5, 3))

# %%
plot_hdi(Exponential())

# %%
plot_hdi(LogNormal(); xlim=(0, 10))

# %%
plot_hdi(Beta(5, 10))

# %%
plot_hdi(Beta(5, 1); legend=:topleft)

# %%
