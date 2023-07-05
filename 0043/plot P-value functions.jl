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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
null = Binomial(30, 1/2)

# %%
pdf(null, 24)

# %%
ccdf(null, 24-1)

# %%
ccdf(Binomial(30, 0.62), 24-1)

# %%
ccdf(Binomial(30, 0.65), 24-1)

# %%
pdf(Binomial(6, 1/6), 3)

# %%
(cdf(Binomial(6, 1/6), 3), ccdf(Binomial(6, 1/6), 3-1))

# %%
safediv(x, y) = ifelse(x == 0, zero(x/y), x/y)

function pvalue_wilson(k, n, p)
    bin = Binomial(n, p)
    μ, σ = mean(bin), std(bin)
    z = safediv(k - μ, σ)
    2ccdf(Normal(), abs(z))
end

function pvalue_clopper_pearson(k, n, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

function pvalue_onetailed_approx(k, n, p; alt = :greater)
    bin = Binomial(n, p)
    μ, σ = mean(bin), std(bin)
    z = safediv(k - μ, σ)
    if alt == :greater
        ccdf(Normal(), z)
    else
        cdf(Normal(), z)
    end
end

function pvalue_onetailed_exact(k, n, p; alt = :greater)
    bin = Binomial(n, p)
    if alt == :greater
        ccdf(bin, k-1)
    else
        cdf(bin, k)
    end
end

function plot_pvalues(; n = 30, k = 21, legend=:topleft)
    @show n, k
    @show pvalue_clopper_pearson(k, n, 0.5)
    @show pvalue_wilson(k, n, 0.5)

    plot(p -> pvalue_clopper_pearson(k, n, p), 0, 1; label="Clopper-Pearson")
    plot!(p -> pvalue_wilson(k, n, p), 0, 1; label="Wilson score", ls=:dot)
    plot!(xtick=0:0.1:1, ytick=0:0.05:1)
    plot!(xguide="parameter p", yguide="P-value")
    title!("binomial data: n = $n,  k = $k"; legend)
end

function plot_onetailed(; n = 30, k = 21, legend=:topleft, alt=:greater)
    @show n, k
    @show alt
    @show pvalue_onetailed_exact(k, n, 0.5; alt)
    @show pvalue_onetailed_approx(k, n, 0.5; alt)

    plot(p -> pvalue_onetailed_exact(k, n, p; alt), 0, 1; label="exact")
    plot!(p -> pvalue_onetailed_approx(k, n, p; alt), 0, 1; label="approx.", ls=:dot)
    plot!(xtick=0:0.1:1, ytick=0:0.05:1)
    plot!(xguide="parameter p", yguide="one-tailed P-value")
    title!("binomial data: n = $n,  k = $k"; legend)
end

# %%
plot_onetailed(; n = 30, k = 21)

# %%
plot_pvalues(; n = 30, k = 21)

# %%
plot_onetailed(; n = 30, k = 9)

# %%
plot_pvalues(; n = 30, k = 9, legend = :topright)

# %%
plot_pvalues(; n = 30, k = 20)

# %%
plot_pvalues(; n = 30, k = 19)

# %%
