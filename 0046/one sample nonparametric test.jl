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

# %%
using Distributions
using HypothesisTests
using QuadGK
using Random
using Roots
using StatsPlots
default(fmt=:png)

cdfsum(dist, x) = quadgk(y -> pdf(dist, y) * cdf(dist, x-y), extrema(dist)...)[1]

empirical_cdf(A, x) = count(≤(x), A) / length(A)

h2(x, y) = (x + y > 0) + (x + y == 0)/2

function nonparametric_t_statistic(X, p=1/2)
    n = length(X)
    phat = 2/(n*(n-1)) * sum(h2(X[i], X[j]) for i in 1:n for j in i+1:n)
    sigmahat2 = 1/(n-1) * sum((mean(h2(X[i], X[j]) for j in 1:n if j != i) - phat)^2 for i in 1:n)
    sehat = 2√(sigmahat2/n)
    t = (phat - p)/sehat
end

function nonparametric_pvalue(X, p=1/2; r=0.6length(X)+1)
    t = nonparametric_t_statistic(X, p)
    pval = 2ccdf(TDist(max(1e-2, n-r)), abs(t))
end

# %%
n = 10
X = rand(Cauchy(), n, 10^5)
T = [nonparametric_t_statistic(x) for x in eachcol(X)]
stephist(T; norm=true, label="T", bin=-6.25:0.5:6.25)
plot!(TDist(0.6n+1); label="TDist($(0.6n+1))")
plot!(Normal(); label="Normal()")
plot!(xlim=(-6, 6))

# %%
n = 15
L = 10^5
X = rand(Cauchy(), n, L)
#X = rand(Normal(), n, L)
pval = [nonparametric_pvalue(x) for x in eachcol(X)]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
#dist = Exponential()
#dist = LogNormal()
dist = InverseGamma(1, 1)
a = find_zero(a -> cdfsum(dist-a, 0) - 0.5, (-10, 10))
b = find_zero(b -> cdfsum(dist, b) - 0.5, (-10, 10))
@show dist a b/2
@show median(dist)
dist = dist - a
n = 20
L = 10^5
X = rand(dist, n, L)
pval = [nonparametric_pvalue(x) for x in eachcol(X)]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
#dist = Exponential()
#dist = LogNormal()
dist = InverseGamma(1, 1)
a = find_zero(a -> cdfsum(dist-a, 0) - 0.5, (-10, 10))
b = find_zero(b -> cdfsum(dist, b) - 0.5, (-10, 10))
@show dist a b/2
@show median(dist)
dist = dist - a
n = 20
L = 10^5
pval = [pvalue(SignedRankTest(rand(dist, n))) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
X = rand(dist, n, L)
T = [nonparametric_t_statistic(x) for x in eachcol(X)]
stephist(T; norm=true, label="T", bin=-6.25:0.5:6.25)
plot!(TDist(0.6n+1); label="TDist($(0.6n+1))")
plot!(Normal(); label="Normal()")
plot!(xlim=(-6, 6))

# %%
plot(dist, -2, 10; label="dist")

# %%
dist1 = LogNormal(log(5), 0.2)
dist2 = LogNormal(log(5), 0.8)
n = 10^7
Δx = rand(dist2, n) - rand(dist1, n)
@show median(Δx)
stephist(Δx; norm=true, label="", xlim=(-10, 40))

# %%
@show median(Δx)
@show m = median(rand(Δx)+rand(Δx) for _ in 1:10^6)
@show m/2;

# %%
n = 100
L = 10^5
pval = [pvalue(SignedRankTest(rand(dist2, n) - rand(dist1, n))) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.05:1)
plot!(size=(400, 400))

# %%
n = 100
L = 10^5
pval = [pvalue(SignedRankTest(rand(dist2-m/2, n) - rand(dist1, n))) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
n = 100
L = 10^5
pval = [nonparametric_pvalue(rand(dist2, n) - rand(dist1, n); ) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.05:1)
plot!(size=(400, 400))

# %%
n = 100
L = 10^5
pval = [nonparametric_pvalue(rand(dist2-m/2, n) - rand(dist1, n); ) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
n = 20
L = 10^5
pval = [pvalue(SignedRankTest(rand(dist2, n) - rand(dist1, n))) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.05:1)
plot!(size=(400, 400))

# %%
n = 20
L = 10^5
pval = [pvalue(SignedRankTest(rand(dist2-m/2, n) - rand(dist1, n))) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
n = 20
L = 10^5
pval = [nonparametric_pvalue(rand(dist2, n) - rand(dist1, n); ) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.05:1)
plot!(size=(400, 400))

# %%
n = 20
L = 10^5
pval = [nonparametric_pvalue(rand(dist2-m/2, n) - rand(dist1, n); ) for _ in 1:L]
plot(x -> empirical_cdf(pval, x), 0, 0.1; label="")
plot!(identity; label="")
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
function nonparametric_phat_sigmahat2(X)
    n = length(X)
    phat = 2/(n*(n-1)) * sum(h2(X[i], X[j]) for i in 1:n for j in i+1:n)
    sigmahat2 = 1/(n-1) * sum((mean(h2(X[i], X[j]) for j in 1:n if j != i) - phat)^2 for i in 1:n)
    phat, sigmahat2
end

function ptilde_sigmatilde2(dist, X)
    ptilde = mean(ccdf(dist, -x) for x in X)
    sigmatilde2 = var(ccdf(dist, -x) for x in X)
    ptilde, sigmatilde2
end

function p_sigma2(dist, n; L=10^5)
    p = 1 - cdfsum(dist, 0)
    Xtmp = zeros(n)
    phat = zeros(L)
    ptilde = zeros(L)
    for i in 1:L
        X = rand!(dist, Xtmp)
        phat[i] = 2/(n*(n-1)) * sum(h2(X[i], X[j]) for i in 1:n for j in i+1:n)
        ptilde[i] = mean(ccdf(dist, -x) for x in X)
    end
    p, mean(phat), n*var(phat)/4, mean(ptilde), n*var(ptilde)
end

dist = Gamma(2, 1) - 2
n = 2^7
L = 10^5
X = rand(dist, n)
[nonparametric_phat_sigmahat2(X), ptilde_sigmatilde2(dist, X), p_sigma2(dist, n; L)]

# %%
