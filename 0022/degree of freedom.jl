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
using LinearAlgebra
using Random
using Distributions
using StatsPlots

# %%
function simlogmlr(dist_true, n, logmaxlik0, logmaxlik1, L=10^5)
    sample = zeros(n)
    logmlr = zeros(L)
    for i in 1:L
        rand!(dist_true, sample)
        logmlr[i] = 2(logmaxlik1(sample) - logmaxlik0(sample))
    end
    logmlr
end

# %%
@time logmlr = simlogmlr(Normal(), 100,
    sample -> loglikelihood(Normal(), sample),
    sample -> loglikelihood(fit_mle(Normal, sample), sample),
)
histogram(logmlr; norm=true, alpha=0.3, bin=0:0.2:15, label="2log(maximum likelihood ratio)")
plot!(Chisq(2), 0, 15; label="Chisq(2)", lw=1.5)
title!("""
    dist_true = Normal(0, 1),  n = 100
    model0 = Normal(0, 1)
    model1 = Normal(μ, σ)
    """; titlefontsize=10)

# %%
@time logmlr = simlogmlr(Normal(), 100,
    sample -> loglikelihood(Normal(mean(sample), 1), sample),
    sample -> loglikelihood(fit_mle(Normal, sample), sample),
)
histogram(logmlr; norm=true, alpha=0.3, bin=0:0.1:6, label="2log(maximum likelihood ratio)")
plot!(Chisq(1), 0.05, 6; label="Chisq(1)", lw=1.5)
title!("""
    dist_true = Normal(0, 1),  n = 100
    model0 = Normal(μ, 1)
    model1 = Normal(μ, σ)
    """; titlefontsize=10)

# %%
@time logmlr = simlogmlr(Normal(), 100,
    sample -> loglikelihood(Normal(0.0, std(sample; corrected=false, mean=0.0)), sample),
    sample -> loglikelihood(fit_mle(Normal, sample), sample),
)
histogram(logmlr; norm=true, alpha=0.3, bin=0:0.1:6, label="2log(maximum likelihood ratio)")
plot!(Chisq(1), 0.05, 6; label="Chisq(1)", lw=1.5)
title!("""
    dist_true = Normal(0, 1),  n = 100
    model0 = Normal(0, σ)
    model1 = Normal(μ, σ)
    """; titlefontsize=10)

# %%
@time logmlr = simlogmlr(Normal(), 100,
    sample -> loglikelihood(Normal(), sample),
    sample -> loglikelihood(Normal(mean(sample), 1), sample),
)
histogram(logmlr; norm=true, alpha=0.3, bin=0:0.1:6, label="2log(maximum likelihood ratio)")
plot!(Chisq(1), 0.05, 6; label="Chisq(1)", lw=1.5)
title!("""
    dist_true = Normal(0, 1),  n = 100
    model0 = Normal(0, 1)
    model1 = Normal(μ, 1)
    """; titlefontsize=10)

# %%
@time logmlr = simlogmlr(Gamma(10, 1), 100,
    sample -> loglikelihood(Gamma(10, 1), sample),
    sample -> loglikelihood(fit_mle(Gamma, sample), sample),
)
histogram(logmlr; norm=true, alpha=0.3, bin=0:0.2:15, label="2log(maximum likelihood ratio)")
plot!(Chisq(2), 0, 15; label="Chisq(2)", lw=1.5)
title!("""
    dist_true = Gamma(10, 1),  n = 100
    model0 = Gamma(10, 1)
    model1 = Gamma(α, θ)
    """; titlefontsize=10)

# %%
n = 5
@time logmlr = simlogmlr(Normal(), n,
    sample -> loglikelihood(Normal(mean(sample), 1.0), sample),
    sample -> sum(x -> logpdf(Normal(x, 1.0), x), sample),
)
histogram(logmlr; norm=true, alpha=0.3, bin=range(0, 4n; length=61), label="2log(maximum likelihood ratio)")
plot!(Chisq(n-1), 0, 4n; label="Chisq(n-1)", lw=1.5)
title!("""
    dist_true = Normal(0, 1),  n = $n
    model0 = MvNormal(fill(μ, n), I)
    model1 = MvNormal(μ⃗, I)
    """; titlefontsize=10)

# %%
n = 5
@time logmlr = simlogmlr(Normal(), n,
    sample -> loglikelihood(Normal(0.0, std(sample; corrected=false, mean=0.0)), sample),
    sample -> loglikelihood(fit_mle(Normal, sample), sample),
)
histogram(logmlr; norm=true, alpha=0.3, bin=0:0.1:6, label="2log(maximum likelihood ratio)")
plot!(Chisq(1), 0.05, 6; label="Chisq(1)", lw=1.5)
title!("""
    dist_true = Normal(0, 1),  n = $n
    model0 = Normal(0, σ)
    model1 = Normal(μ, σ)
    """; titlefontsize=10)

# %%
n = 5
@time logmlr = simlogmlr(Normal(), n,
    sample -> loglikelihood(Normal(0.0, std(sample; corrected=false, mean=0.0)), sample),
    sample -> loglikelihood(fit_mle(Normal, sample), sample),
)
fstat = @. (n - 1)*(exp(logmlr/n) - 1)

P = histogram(fstat; norm=true, alpha=0.3, bin=range(0, 6; length=61), label="(t-statistic)²")
plot!(FDist(1, n-1), 0.05, 6; label="FDist(1, n-1)", lw=1.5)
plot!(Chisq(1), 0.05, 6; label="Chisq(1)", lw=1.5, ls=:dot)
title!("dist_true = Normal(0, 1),  n = $n"; titlefontsize=10)

Q = histogram(fstat; norm=true, alpha=0.3, bin=range(0, 30; length=61), label="(t-statistic)²")
plot!(FDist(1, n-1), 4, 30; label="FDist(1, n-1)", lw=1.5)
plot!(Chisq(1), 4, 30; label="Chisq(1)", lw=1.5, ls=:dash)
title!("dist_true = Normal(0, 1),  n = $n"; titlefontsize=10)
plot!(; xlim=(4, 30), ylim=(0, 0.034))

plot(P, Q; size=(800, 300))

# %%
function simfstat(dist_true, n, L=10^5)
    sample = zeros(n)
    fstat = zeros(L)
    μ₀ = mean(dist_true)
    σ₀² = var(dist_true)
    for i in 1:L
        rand!(dist_true, sample)
        fstat[i] = (mean(sample) - μ₀)^2/(var(sample)/n)
    end
    fstat
end

n = 5
@time fstat = simfstat(Normal(2, 1), n)

P = histogram(fstat; norm=true, alpha=0.3, bin=range(0, 6; length=61), label="(t-statistic)²")
plot!(FDist(1, n-1), 0.05, 6; label="FDist(1, n-1)", lw=1.5)
plot!(Chisq(1), 0.05, 6; label="Chisq(1)", lw=1.5, ls=:dot)
title!("dist_true = Normal(0, 1),  n = $n"; titlefontsize=10)

Q = histogram(fstat; norm=true, alpha=0.3, bin=range(0, 30; length=61), label="(t-statistic)²")
plot!(FDist(1, n-1), 4, 30; label="FDist(1, n-1)", lw=1.5)
plot!(Chisq(1), 4, 30; label="Chisq(1)", lw=1.5, ls=:dash)
title!("dist_true = Normal(0, 1),  n = $n"; titlefontsize=10)
plot!(; xlim=(4, 30), ylim=(0, 0.034))

plot(P, Q; size=(800, 300))

# %%
