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
# * 2022-08-06
# * License: https://opensource.org/licenses/MIT

# %%
using Distributions
using Roots
using StatsPlots
default(fmt=:png, size=(400, 250))

# %%
"""
    invpdf(dist, y; alg = Order0(), modefunc = mymode)

`dist` は最頻値を持つ単峰型の連続分布であると仮定する.

この函数はその分布の最頻値の左側と右側での pdf の逆函数の `y` における値達のタプルを返す.
"""
function invpdf(dist, y; alg = Order0(), modefunc = mymode)
    m = modefunc(dist)
    cdfm = cdf(dist, m)
    a0 = quantile(dist, cdfm/2)
    b0 = quantile(dist, 1 - (1 - cdfm)/2)
    f(x) = logpdf(dist, x) - log(y)
    a = m == minimum(dist) ? minimum(dist) : find_zero(f, a0, alg)
    b = m == maximum(dist) ? maximum(dist) : find_zero(f, b0, alg)
    a, b
end

"""
    mymode(x)

適当に修正された `mode(x)`.
"""
mymode(x) = mode(x)

function mymode(beta::Beta)
    α, β = params(beta)
    α == 1 && β == 1 && return median(beta)
    α == 1 && β > 1 && return mininmum(beta)
    α > 1 && β == 1 && return maximum(beta)
    mode(beta)
end

@doc invpdf

# %%
@doc mymode

# %%
dist = Normal()
invpdf(dist, pdf(dist, mode(dist))/2)

# %%
dist = Gamma(2, 3)
invpdf(dist, pdf(dist, mode(dist))/2)

# %%
dist = Gamma(2, 3)
invpdf(dist, 0.05)

# %%
dist = Exponential()
invpdf(dist, 0.1)

# %%
dist = Beta(5, 1)
invpdf(dist, 0.1)

# %%
"""
    cdfinvpdf(dist, y; alg = Order0())

`dist` は `mode(dist)` を持つ単峰型の連続分布であると仮定する.

この函数はその分布の左側と右側での pdf の逆函数の `y` における値達のあいだの区間の確率の値を返す.

"""
function cdfinvpdf(dist, y; alg = Order0())
    a, b = invpdf(dist, y; alg)
    cdf(dist, b) - cdf(dist, a)
end

@doc cdfinvpdf

# %%
dist = Gamma(2, 3)
m = mode(dist)
pdfm = pdf(dist, m)
plot(y -> cdfinvpdf(dist, y), eps(), pdfm)

# %%
"""
    hdi(dist, α = 0.05; alg = Order0(), modefunc = mymode)

`dist` は `mode(dist)` を持つ単峰型の連続分布であると仮定する.

この函数はその分布の100(1-α)% HDI (highest density interval)を返す.
"""
function hdi(dist, α = 0.05; alg = Order0(), modefunc = mymode)
    pdfm = pdf(dist, modefunc(dist))
    y = find_zero(pdfm/2, alg) do y
        cdfinvpdf(dist, y; alg) - (1 - α)
    end
    invpdf(dist, y; alg)
end

@doc hdi

# %%
function plot_hdi(dist, α = 0.05; alg=Order0(), kwargs...)
    @show α
    @show a, b = hdi(dist, α; alg)
    plot(dist; label="pdf")
    vline!([a, b]; label="hdi")
    plot!(; kwargs...)
end

# %%
plot_hdi(Normal())

# %%
plot_hdi(Gamma(10, 3))

# %%
plot_hdi(Exponential())

# %%
plot_hdi(Gamma(10, 1))

# %%
plot_hdi(Beta(5, 10))

# %%
plot_hdi(Beta(5, 1); legend=:topleft)

# %%
