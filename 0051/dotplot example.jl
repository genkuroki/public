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
#     display_name: Julia 1.11.1
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# https://qiita.com/JMP_Japan/items/4a924c361206e2b862cd

# %%
using Distributions
using StatsPlots
default(fmt=:png)

function pvalue_onesample_t_test(n, x̄, s², μ=0.0)
    z = (x̄ - μ) / √(s² / n)
    2ccdf(TDist(n - 1), abs(z))
end

function pvalue_onesample_t_test(x, μ=0.0)
    n, x̄, s² = length(x), mean(x), var(x)
    pvalue_onesample_t_test(n, x̄, s², μ)
end

function confint_onesample_t_test(n, x̄, s², α=0.05)
    c = cquantile(TDist(n - 1), α / 2)
    [x̄ - c * √(s² / n), x̄ + c * √(s² / n)]
end

function confint_onesample_t_test(x, α=0.05)
    n, x̄, s² = length(x), mean(x), var(x)
    confint_onesample_t_test(n, x̄, s², α)
end

# %%
data = Any[
    1   163.5
    2   177.8
    3   183.0
    4   172.7
    5   168.4
    6   175.5
    7   160.3
    8   180.2
    9   165.8
    10  174.6
    11  162.9
    12  177.4
    13  180.5
    14  169.8
    15  176.5
]

x = oftype.(data[1, 2], data[:, 2])
@show x
@show sort(x)
@show quantile.((x,), (0.00, 0.25, 0.50, 0.75, 1.00))

μ = 168.3
n, x̄, s = length(x), mean(x), std(x)
@show n x̄ s
@show pvalue_onesample_t_test(x, μ)
@show confint_onesample_t_test(x)

dotplot(["x";;], [x;;]; label="", msc=:auto, ms=3)
plot!(xlim=(-0.5, 1.5), ylim=(158, 184))
plot!(ytick=150:200, tickfontsize=6)
plot!(size=(100, 300))

# %%
r = 14
X = rand(Normal(x̄, s), n, r)
xs = "x" .* string.((1:r)')
dotplot(["x" xs], [x X]; label="", msc=:auto, ms=3)
plot!(ytick=150:2:200, tickfontsize=6)
plot!(size=(1000, 300))

# %%
