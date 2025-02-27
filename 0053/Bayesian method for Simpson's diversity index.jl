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
#     display_name: Julia 1.11.3
#     language: julia
#     name: julia-1.11
# ---

# %%
using Distributions
using KernelDensity
using LinearAlgebra
using Optim
using QuadGK
using Random
using Roots
using StatsPlots
default(fmt=:png)
using BenchmarkTools: @btime, @benchmark

simpson_diversity_index(p) = 1 - dot(p, p)

function sidhat_naive(K)
    n = sum(K)
    1 - sum(k*(k-1)/(n*(n-1)) for k in K)
end

function sidhat(K)
    n = sum(K)
    1 - (dot(K, K) - n)/(n*(n-1))
end

sid_bayes_naive(P, n) = 1 - sum(p^2 - 2p*(1-p)/(n-1) for p in P)
sid_bayes(P, n) = 1  - ((1 + 2/(n-1))*dot(P, P) - 2/(n-1))

function sehat_sdi(x)
    n = sum(x)
    √max(0, (4/n) * (sum(x->(x/n)^3, x) - sum(x -> (x/n)^2, x)^2))
end

function confint_sdi(x; α=0.05)
    D̂ = sidhat(x)
    sehat = sehat_sdi(x)
    w = cquantile(Normal(), α/2)
    [D̂ - w*sehat, D̂ + w*sehat]
end

r(x) = round(x; sigdigits=3)

function plot_posterior(; data=[5, 5, 5, 5], κ=0, α=0.05, L=10^6, kwargs...)
    m = length(data)
    D̂ = sidhat(data)
    confint_D = confint_sdi(data)
    prior_param = fill(κ, m)
    posterior = Dirichlet(data .+ prior_param)
    n = sum(params(posterior)[1])
    Ptmp = zeros(m)
    D = [sid_bayes(rand!(posterior, Ptmp), n) for _ in 1:L]
    credint_D = quantile.((D,), [α/2, 1-α/2])
    mean_D = mean(D)
    median_D = median(D)
    ik_D = InterpKDE(kde(D))
    mode_D = optimize(d -> -pdf(ik_D, d), minimum(D), maximum(D)).minimizer
    
    @show data
    @show prior_param

    density(D; label="posterior")
    plot!(credint_D, zeros(2); label="95% cred. int. = $(r.(credint_D))", lw=3)
    scatter!([mean_D], [0.0]; label="mean = $(r(mean_D))")
    #scatter!([median_D], [0.0]; label="median = $(r(median_D))")
    #scatter!([mode_D], [pdf(ik_D, mode_D)]; label="mode = $(r(mode_D))")
    plot!(confint_D, fill(-0.05pdf(ik_D, mode_D), 2); ls=:dash, 
        label="95% conf.int = $(r.(confint_D))", lw=3)
    #scatter!([D̂], [-0.05pdf(ik_D, mode_D)]; label="sample SDI = $(r(D̂))")
    vline!([D̂]; label="sample SDI = $(r(D̂))")
    plot!(xguide="Simpson's diversity index", yguide="posterior density")
    plot!(; kwargs...)
    #title!("data=$data, prior=$prior_param")
end

# %%
K = 10rand(100)
@btime sidhat($K)
@btime sidhat_naive($K)
@show sidhat(K) ≈ sidhat_naive(K)

n = sum(K)
P = K / n
@btime sid_bayes($P, n)
@btime sid_bayes_naive($P, n)
@show sid_bayes(P, n) ≈ sid_bayes_naive(P, n)
;

# %%
n = 10
p = rand(Exponential(), n)
p ./= sum(p)
mult = Multinomial(n, p)
Ktmp = rand(mult)
@time D̂ = [sidhat(rand!(mult, Ktmp)) for _ in 1:10^5]
@show simpson_diversity_index(p)
@show mean(D̂)
@time SEhat = [sehat_sdi(rand!(mult, Ktmp)) for _ in 1:10^5]
@show std(D̂)
@show mean(SEhat)
@show mean(SEhat) * n/(n-1)
;

# %%
plot_posterior(data = [5, 5, 5, 5])

# %%
plot_posterior(data = [1, 6])

# %% [markdown]
# https://journals.asm.org/doi/10.1128/jcm.39.11.4190-4192.2001

# %%
table1 = [9; 8; 7; 6; 5; 5; fill(4, 3); fill(3, 4); fill(2, 9); fill(1, 35)]
plot_posterior(data=table1, legend=:topleft)

# %%
table2 = [37; 5; fill(4, 2); fill(3, 4); fill(2, 8); fill(1, 39)]
plot_posterior(data=table2)

# %%
table3 = [30; 13; 9; 8; fill(7, 3); fill(6, 2); 5; fill(2, 3); fill(1, 13)]
plot_posterior(data=table3)

# %%
function pvalues_bayes(data, sid₀; κ=eps(), L=10^4)
    m = length(data)
    n = sum(data)
    prior_param = fill(κ, m)
    posterior = Dirichlet(data .+ prior_param)
    Ptmp = rand(posterior)
    c = zeros(Int, length(sid₀))
    d = zeros(Int, length(sid₀))
    for i in 1:L
        P = rand!(posterior, Ptmp)
        sid = sid_bayes(P, n)
        @. c += sid ≤ sid₀
        @. d += sid ≥ sid₀
    end
    @. min(1, 2c/L, 2d/L)
end

function pvalue_bayes(data, sid₀; κ=eps(), L=10^4)
    m = length(data)
    n = sum(data)
    prior_param = fill(κ, m)
    posterior = Dirichlet(data .+ prior_param)
    Ptmp = rand(posterior)
    c = 0
    d = 0
    for i in 1:L
        P = rand!(posterior, Ptmp)
        sid = sid_bayes(P, n)
        c += sid ≤ sid₀
        d += sid ≥ sid₀
    end
    min(1, 2c/L, 2d/L)
end

function sim_pval(dist::Multinomial; niters=10^4, κ=eps(), L=10^4)
    μ = mean(dist)
    sdi₀ = sidhat(μ)
    pval = zeros(niters)
    Ktmp = [rand(dist) for _ in 1:Threads.nthreads()]
    Threads.@threads :static for i in 1:niters
        tid = Threads.threadid()
        K = rand!(dist, Ktmp[tid])
        pval[i] = pvalue_bayes(K, sdi₀; κ, L)
    end
    pval
end

_ecdf(A, x) = count(≤(x), A) / length(A)

# %%
data = table1

@show D̂ = sidhat(data)
sid₀ = range(D̂-0.02, D̂+0.02, 100)
@show length(sid₀)
@time pval = pvalues_bayes(data, sid₀; κ=eps(), L=10^6)
plot(sid₀, pval)
vline!([D̂])

# %%
data = table1
dist = Multinomial(sum(data), data/sum(data))
@time pval = sim_pval(dist; κ=0.23, niters=10^4)
plot(α -> _ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot)
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
data = table2
dist = Multinomial(sum(data), data/sum(data))
@time pval = sim_pval(dist; κ=0.05, niters=10^4)
plot(α -> _ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot)
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
data = table3
dist = Multinomial(sum(data), data/sum(data))
@time pval = sim_pval(dist; κ=0.1, niters=10^4)
plot(α -> _ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot)
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
data = [10, 10, 10, 10]
dist = Multinomial(sum(data), data/sum(data))
@time pval = sim_pval(dist; κ=0.1, niters=10^4)
plot(α -> _ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot)
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
data = [10, 6, 2, 1, 1]
dist = Multinomial(sum(data), data/sum(data))
@time pval = sim_pval(dist; κ=0.1, niters=10^4)
plot(α -> _ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot)
plot!(xtick=0:0.01:1, ytick=0:0.01:1)
plot!(size=(400, 400))

# %%
