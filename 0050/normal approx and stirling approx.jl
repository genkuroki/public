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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png, legendfontsize=9, size=(500, 350))

n, p = 10, 0.3
@show n*p
bin = Binomial(n, p)
normal = Normal(mean(bin), std(bin))
bar(bin; alpha=0.3, label="Binomial(n=$n, p=$p)")
plot!(normal; lw=2.5, label="normal approx.")
plot!(xtick=0:n)
title!("normal approx. of binomial distribution") |> display

using Plots
default(fmt=:png, legendfontsize=10, size=(500, 350))

stirling_aporox(n) = n^n * exp(-n) * √(2π*n)
nmax = 6
scatter(1:nmax, n -> factorial(n); label="n!")
plot!(stirling_aporox, 0.5, nmax+0.5; label="Stirling approx.")
ytick = 2 .^ (0:11)
plot!(yscale=:log10, ytick=(ytick, string.(ytick)))

# %%

# %%
