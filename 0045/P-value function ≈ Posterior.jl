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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using Optim
using Roots
using StatsPlots
default(fmt=:png, titlefontsize=10)

function pvalue_wilson(k, n, p)
    p̂ = k/n
    SE = √(p*(1-p)/n)
    2ccdf(Normal(), abs(p̂ - p)/SE)
end

function confint_wilson(k, n; α = 0.05)
    p̂ = k/n
    z = quantile(Normal(), 1-α/2)
    a, b, c = 1+z^2/n, p̂+z^2/(2n), p̂^2
    # ap² - 2bp + c = 0 を解く.
    sqrtD = √(b^2 - a*c)
    p_L = (b - sqrtD)/a
    p_U = (b + sqrtD)/a
    [p_L, p_U]
end

function pvalue_hdi(dist::ContinuousUnivariateDistribution, x₀; xlim = extrema(dist))
    p₀ = pdf(dist, x₀)
    m = mode(dist)
    f(x) = pdf(dist, x) - p₀
    if x₀ == m
        1.0
    elseif x₀ > m
        x₁ = find_zero(f, (xlim[begin], m))
        cdf(dist, x₁) + ccdf(dist, x₀)
    else
        x₁ = find_zero(f, (m, xlim[end]))
        cdf(dist, x₀) + ccdf(dist, x₁)
    end
end

function pvalue_bayes_hdi(k, n, p; a=1, b=1)
    posterior = Beta(k+a, n-k+b)
    if k+a ≤ 1
        return ccdf(posterior, p)
    elseif n-k+b ≤ 1
        return cdf(posterior, p)
    end
    pvalue_hdi(posterior, p)
end

function hdi(dist::ContinuousUnivariateDistribution, α = 0.05; alg = Brent())
    f(t) = quantile(dist, t + (1 - α)) - quantile(dist, t)
    o = optimize(f, 0, α, alg)
    p = o.minimizer
    quantile.(dist, [p, p + (1 - α)])
end

function confint_bayes_hdi(k, n; α=0.05, a=1, b=1)
    posterior = Beta(k+a, n-k+b)
    hdi(posterior, α)
end

function plot_pvalue_functions_and_posterior(k, n; a=1, b=1)
    ps = range(0, 1, 1001)
    
    P = plot()
    plot!(ps, p -> pvalue_wilson(k, n, p); label="Wilson")
    plot!(ps, p -> pvalue_bayes_hdi(k, n, p; a, b); label="Bayesian", ls=:dash)
    plot!(xtick=0:0.05:1, ytick=0:0.1:1)
    title!("P-value functions for n=$n, k=$k")
    
    posterior = Beta(k+a, n-k+b)
    A = pdf(posterior, k+a ≤ 1 ? 0.0 : n-k+b ≤ 1 ? 1.0 : mode(posterior))
    Q = plot()
    plot!(ps, p -> pdf(posterior, p)/A; label="", c=2)
    plot!(xtick=0:0.05:1, ytick=0:0.1:1)
    title!("posterior pdf normaluzed by max=1 for n=$n, k=$k")

    plot(P, Q; size=(600, 400), layout=(2, 1))
    plot!(titlefontsize=10, tickfontsize=6)
end

# %%
plot_pvalue_functions_and_posterior(3, 10)

# %%
plot_pvalue_functions_and_posterior(6, 20)

# %%
N = 2000
_ns = round.(Int, exp.(range(log(5), log(N), 400)))
ns = [fill(_ns[1], 20); _ns; fill(_ns[end], 20)]

p0 = 0.3
K = cumsum(rand(Bernoulli(p0), N))
@gif for n in ns
    plot_pvalue_functions_and_posterior(K[n], n)
end

# %%
