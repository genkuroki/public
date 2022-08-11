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

# %%
using Distributions
using StatsPlots
default(fmt=:png, size=(500, 280), titlefontsize=10, tickfontsize=6)
using SymPy

function dist_of_sum_of_dice_rolls(n; r = 6)
    @vars z
    f = sum(z^i for i in 1:r)
    fr = expand(f^n)
    xs = n:n*r
    ys = N.(sympy.Poly(fr).coeffs())
    ps = ys / sum(ys)
    f, fr, xs, ps
end

function normal_approx_of_dist_of_sum_of_dice_rolls(n; r = 6)
    μ = (1 + r)/2
    σ² = mean((i - μ)^2 for i in 1:r)
    Normal(n*μ, √(n*σ²))
end

function plot_dist_of_sum_of_dice_rolls(n; r = 6, xtick=nothing)
    f, fr, xs, ps = dist_of_sum_of_dice_rolls(n; r)
    normal = normal_approx_of_dist_of_sum_of_dice_rolls(n; r)
    display(f)
    display(fr)
    @show normal
    isnothing(xtick) && (xtick = xs)
    bar(xs, ps; label="", title="n = $n", alpha=0.4, xtick)
    plot!(normal; label="", lw=1.5)
end

# %%
plot_dist_of_sum_of_dice_rolls(1)

# %%
plot_dist_of_sum_of_dice_rolls(2)

# %%
plot_dist_of_sum_of_dice_rolls(3)

# %%
plot_dist_of_sum_of_dice_rolls(4)

# %%
plot_dist_of_sum_of_dice_rolls(5)

# %%
plot_dist_of_sum_of_dice_rolls(6)

# %%
plot_dist_of_sum_of_dice_rolls(7; xtick=7:2:7*6)

# %%
plot_dist_of_sum_of_dice_rolls(8; xtick=7:2:7*6)

# %%
plot_dist_of_sum_of_dice_rolls(9; xtick=9:2:9*6)

# %%
plot_dist_of_sum_of_dice_rolls(10; xtick=10:2:10*6)

# %%
