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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# # 二項分布モデルでのP値函数のグラフ

# %% [markdown]
# ## WilsonのP値函数の場合

# %%
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6, guidefontsize=9)

function pvalue_bin_wilson(k, p; n=20)
    z = (k - n*p)/√(n*p*(1 - p))
    2ccdf(Normal(), abs(z))
end

# %%
n = 20
α = 0.05
x = range(-0.5, n+0.5, 500)
p = range(0, 1, 500)
heatmap(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n) ≥ α;
    clim=(0, 1), colorbar=false)
plot!(xguide="data k", yguide="parameter p")
plot!(xtick=0:n, ytick=0:0.1:1)
title!("Wilson's confidence intervals for n = $n, α = $α")
plot!(size=(580, 500))

# %%
n = 20
x = range(-0.5, n+0.5, 500)
p = range(0, 1, 500)
heatmap(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); clim=(0, 1))
plot!(xguide="data k", yguide="parameter p")
plot!(xtick=0:n, ytick=0:0.1:1)
title!("Wilson's P-value function for n = $n")
plot!(size=(640, 500))

# %%
n, k = 20, 6
p = range(0, 1, 500)
plot(p, p -> pvalue_bin_wilson(k, p; n); label="")
plot!(xguide="parameter p", yguide="P-value")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
title!("Wilson's P-value function for n=$n, k=$k")

# %%
n = 20
x = range(-0.5, n+0.5, 211)
p = range(0, 1, 201)
@time anim = @animate for t in 0:3:359
    surface(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); colorbar=false, camera=(t, 60))
    plot!(xguide="k", yguide="p", zguide="P-value")
    plot!(size=(640, 600), zlim=(0, 1))
end
gif(anim, "Wilson's P-value function n=$n (60).gif")

# %%
n = 20
x = range(-0.5, n+0.5, 211)
p = range(0, 1, 201)
@time anim = @animate for t in 0:3:359
    surface(x, p, (x, p)->pvalue_bin_wilson(round(x), p; n); colorbar=false, camera=(t, 45))
    plot!(xguide="k", yguide="p", zguide="P-value")
    plot!(size=(640, 600), zlim=(0, 1))
end
gif(anim, "Wilson's P-value function n=$n (45).gif")

# %% [markdown]
# ## WilsonのP値函数とClopper-PeasonのP値函数とベイズ版P値函数の比較

# %%
using Distributions
using StatsPlots
default(fmt=:png,
    titlefontsize=10, tickfontsize=6, guidefontsize=9, legendfontsize=9)
safediv(x, y) = x==0 ? x : y==0 ? oftype(x, Inf) : x/y

# WilsonのP値函数
function pvalue_bin_wilson(k, p; n=20)
    z = safediv(k - n*p, √(n*p*(1 - p)))
    2ccdf(Normal(), abs(z))
end

# Clopper-PearsonのP値函数
function pvalue_bin_cp(k, p; n=20)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

# 事前分布Beta(γ, δ)によるBayes版のP値函数(equal-tailed interval版)
function pvalue_bin_bayesian_eti(k, p; n=20, γ=1, δ=1)
    beta = Beta(k+γ, n-k+δ)
    min(1, 2cdf(beta, p), 2ccdf(beta, p))
end

function plot_pvalue_functions(; n=20, k=6, γ=1, δ=1, kwargs...)
    p = range(0, 1, 500)
    plot(p, p -> pvalue_bin_wilson(k, p; n); label="Wilson")
    plot!(p, p -> pvalue_bin_cp(k, p; n); label="Clopper-Pearson", ls=:dash)
    plot!(p, p -> pvalue_bin_bayesian_eti(k, p; n, γ, δ); label="Bayesan", ls=:dashdot)
    plot!(xguide="parameter p", yguide="P-value")
    plot!(xtick=0:0.1:1, ytick=0:0.1:1)
    title!("P-value functions for n=$n, k=$k, prior=Beta($γ, $δ)")
    plot!(; kwargs...)
end

# %%
plot_pvalue_functions(; n=20, k=6, γ=1, δ=1)

# %%
PP = []
for k in 0:11
    P = plot_pvalue_functions(; n=20, k, title="n=$n, k=$k", legend=false)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(4, 3))
plot!(legend=false)

# %%
plot_pvalue_functions(; n=20, k=6, γ=1//3, δ=1//3)

# %%
PP = []
for k in 0:11
    P = plot_pvalue_functions(; n=20, k, γ=1//3, δ=1//3, title="n=$n, k=$k", legend=false)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(4, 3))
plot!(legend=false)

# %%
