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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using Distributions
using HypothesisTests
using StatsPlots
default(size=(500, 300))

function pval_t(null, X)
    μ = mean(null)
    n = length(X)
    T = (mean(X) - μ)/√(var(X)/n)
    2ccdf(TDist(n-1), abs(T))
end

function pval_wilcoxon(null, X)
    m = quantile(null, 0.5) # median of null
    W = HypothesisTests.SignedRankTest(X .- m)
    pvalue(W)
end

function plot_pvals(;
        pval = pval_wilcoxon,
        null = Gamma(4, 1),
        n = 100,
        L = 10^6
    )
    PVal = similar(zeros(), L)
    for i in 1:L
        X = rand(null, n)
        PVal[i] = pval(null, X)
    end
    realalpha = mean(PVal .< 0.05)
    histogram(PVal; norm=true, alpha=0.3, bin=0.001:0.025:1.001, label="")
    plot!(; xtick=0:0.05:1, xrotation=90, bottom_margin=3Plots.mm)
    title!("""
        pval = $pval,  n = $n
        null = $null
        P(p-value < 0.05) = $realalpha,  niters = $L""",
        titlefontsize=8)
end

function plot_both(;
        null = Gamma(4, 1),
        n = 100,
        L = 10^6
    )
    P = plot_pvals(; pval = pval_t, null, n, L)
    Q = plot_pvals(; pval = pval_wilcoxon, null, n, L)
    plot(P, Q; size=(500, 600), layout=(2, 1), tickfontsize=8)
end

# %%
plot(Gamma(4, 1); label="Gamma(4, 1)")

# %%
plot_both(; null = Gamma(4, 1), n = 30)

# %%
plot_both(; null = Gamma(4, 1), n = 100)

# %%
plot_both(null = Gamma(4, 1), n = 300, L=10^5)

# %%
plot_both(; null = Gamma(4, 1), n = 1000, L=10^5)

# %%
plot(Gamma(10, 1); label="Gamma(10, 1)")

# %%
plot_both(; null = Gamma(10, 1))

# %%
plot(Gamma(20, 1); label="Gamma(20, 1)")

# %%
plot_both(; null = Gamma(20, 1))

# %%
plot(Gamma(50, 1), 20, 80; label="Gamma(50, 1)")

# %%
plot_both(; null = Gamma(50, 1), n = 30)

# %%
plot_both(; null = Gamma(50, 1), n = 100)

# %%
plot_both(; null=Gamma(50, 1), n = 300, L=10^5)

# %%
plot_both(; null=Gamma(50, 1), n = 1000, L=10^5)

# %%
plot(Uniform(-1, 1), -3, 3; label="Uniform(-1, 1)")

# %%
plot_both(; null = Uniform(-1, 1))

# %%
null = TDist(3)
plot(null, -8, 8; label="TDist(3)")
plot!(Normal(mean(null), std(null)), -8, 8; label="normal approx.", ls=:dash)

# %%
plot_both(; null)

# %%
null = MixtureModel([Normal(), Normal(10, 1)], [0.95, 0.05])
plot(null, -5, 15; label="", title="$null", titlefontsize=8) |> display
plot_both(; null, n = 30)

# %%
plot_both(; null, n = 100)

# %%
plot_both(; null, n = 300, L = 10^5)

# %%
plot_both(; null, n = 1000, L = 10^5)

# %%
null = Exponential()
plot(null; label="", title="$null", titlefontsize=8) |> display
plot_both(; null, n = 30)

# %%
plot_both(; null, n = 100)

# %%
