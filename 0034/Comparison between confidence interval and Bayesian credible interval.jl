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
default(fmt=:png, titlefontsize=10, size=(400, 250))

# %%
function pvalue_clopper_pearson(n, k, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

function confint_clopper_pearson(n, k; α=0.05)
    p_L = k > 0 ? quantile(Beta(k, n-k+1), α/2) : zero(α)
    p_U = k < n ? quantile(Beta(k+1, n-k), 1-α/2) : one(α)
    [p_L, p_U]
end

function pvalue_wilson(n, k, p)
    p̂ = k/n
    SE = √(p*(1-p)/n)
    2ccdf(Normal(), abs(p̂ - p)/SE)
end

function confint_wilson(n, k; α=0.05)
    p̂ = k/n
    z = quantile(Normal(), 1-α/2)
    a, b, c = 1+z^2/n, p̂+z^2/(2n), p̂^2
    # ap² - 2bp + c = 0 を解く.
    sqrtD = √(b^2 - a*c)
    p_L = (b - sqrtD)/a
    p_U = (b + sqrtD)/a
    [p_L, p_U]
end

function pvalue_bayesian(n, k, p; prior=(1/2, 1/2))
    beta = Beta((prior .+ (k, n-k))...)
    min(1, 2cdf(beta, p), 2ccdf(beta, p))
end

function credint_bayesian(n, k; α=0.05, prior=(1/2, 1/2))
    beta = Beta((prior .+ (k, n-k))...)
    quantile.(beta, [α/2, 1-α/2])
end

function show_cis(n, k; α=0.05, prior=(1/2, 1/2),
        s=Bool[1,1,0], legend=:topleft)
    ci_ba = credint_bayesian(n, k; α, prior)
    ci_wi = confint_wilson(n, k; α)
    ci_cp = confint_clopper_pearson(n, k; α)
    beta0 = Beta(prior...)
    beta = Beta((prior .+ (k, n-k))...)
    prob_ci_ba = cdf(beta, ci_ba[end]) - cdf(beta, ci_ba[begin])
    prob_ci_wi = cdf(beta, ci_wi[end]) - cdf(beta, ci_wi[begin])
    prob_ci_cp = cdf(beta, ci_cp[end]) - cdf(beta, ci_cp[begin])
    if s[1]
        @show ci_ba
        @show prob_ci_ba
    end
    if s[2]
        @show ci_wi
        @show prob_ci_wi
    end
    if s[3]
        @show ci_cp
        @show prob_ci_cp
    end
    
    P1 = plot(; legend)
    plot!(beta; label="posterior", c=1)
    plot!(beta0; label="prior", ls=:dot, c=:black)
    title!("prior and posterior of data (n=$n, k=$k)")
    plot!(tickfontsize=6, xtick=0:0.1:1)
    
    P2 = plot(; legend)
    if s[1]
        plot!(p -> pvalue_bayesian(n, k, p; prior); label="Bayes", c=1)
        plot!(ci_ba, fill(α+0.005, 2); label="", c=1)
    end
    if s[2]
        plot!(p -> pvalue_wilson(n, k, p); label="Wilson", c=2, ls=:dash)
        plot!(ci_wi, fill(α, 2); label="", c=2, ls=:dash)
    end
    if s[3]
        plot(p -> pvalue_clopper_pearson(n, k, p), 0, 1; label="CP", c=3, ls=:dashdot)
        plot!(ci_cp, fill(α-0.005, 2); label="", c=3, ls=:dashdot)
    end
    title!("P-value functions of data (n=$n, k=$k)")
    plot!(tickfontsize=6, xtick=0:0.1:1, ytick=0:0.05:1)

    plot(P1, P2; size=(800, 250))
end

# %%
prior = (7, 3)
beta0 = Beta(prior...)
plot(beta0; label="", title="prior: $beta0", ls=:dot, c=:black)

# %%
show_cis(2, 1; prior=(7, 3))

# %%
show_cis(4, 2; prior=(7, 3))

# %%
show_cis(10, 5; prior=(7, 3))

# %%
show_cis(20, 10; prior=(7, 3))

# %%
show_cis(40, 20; prior=(7, 3))

# %%
show_cis(80, 40; prior=(7, 3))

# %%
show_cis(160, 80; prior=(7, 3))

# %%
show_cis(320, 160; prior=(7, 3))

# %%
function prob_ci_wilson(n, k; α=0.05, prior=(1/2, 1/2))
    beta = Beta((prior .+ (k, n-k))...)
    ci = confint_wilson(n, k; α)
    cdf(beta, ci[end]) - cdf(beta, ci[begin])
end

function plot_prob_ci(; α=0.05, p_true=0.5, prior=(1/2, 1/2),
        nmin=2, nmax=1000, kwargs...)
    plot(nmin:2:nmax, n -> prob_ci_wilson(n, n*p_true; α, prior); label="")
    title!("prior = $prior,  data = (n, $(p_true)n),  α = $α\n\
        probability of p ∈ Wilson's CI with p ∼ posterior")
    hline!([1-α]; label="", ls=:dot)
    plot!(xguide="n")
    plot!(topmargin=4Plots.mm)
    plot!(; kwargs...)
end

# %%
plot_prob_ci(prior=(7, 3), nmax=1000, ylim=(0.93, 0.97))

# %%
plot_prob_ci(prior=(1, 1), nmax=1000, ylim=(0.93, 0.97))

# %%
show_cis(10, 5; prior=(7, 3))

# %%
show_cis(10, 5; prior=(1, 1))

# %%
show_cis(10, 7; prior=(7, 3))

# %%
show_cis(10, 7; prior=(1, 1))

# %%
plot_prob_ci(prior=(7, 3), nmax=1000, ylim=(0.93, 0.97), p_true=0.7)

# %%
plot_prob_ci(prior=(1, 1), nmax=1000, ylim=(0.93, 0.97), p_true=0.7)

# %%
show_cis(10, 3; prior=(7, 3), legend=:topright)

# %%
show_cis(10, 3; prior=(1, 1), legend=:topright)

# %%
plot_prob_ci(prior=(7, 3), nmax=1000, ylim=(0.8, 0.97), p_true=0.3)

# %%
plot_prob_ci(prior=(1, 1), nmax=1000, ylim=(0.93, 0.97), p_true=0.3)

# %%
