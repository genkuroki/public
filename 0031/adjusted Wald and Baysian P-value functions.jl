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
using Base.Threads
using Distributions
using Random
using StatsPlots
default(fmt=:png, titlefontsize=10)

# %%
safediv(x, y) = x == 0 ? x : x/y

function pvalue_clopper_pearson(n, k, p)
    k == 0 && return min(1, 2ccdf(Beta(k+1, n-k), p))
    k == n && return min(1, 2cdf(Beta(k, n-k+1), p))
    min(1, 2cdf(Beta(k, n-k+1), p), 2ccdf(Beta(k+1, n-k), p))
end

function pvalue_wilson_score(n, k, p)
    2ccdf(Normal(), safediv(abs(k - n*p), √(n*p*(1-p))))
end

function pvalue_adjusted_wald(n, k, p; z = 2)
    p̂ = (k + z^2/2) / (n + z^2)
    2ccdf(Normal(), safediv(abs(p̂ - p), √(p̂*(1-p̂)/n)))
end

function pvalue_bayesian(n, k, p; a = 1, b = 1)
    min(1, 2cdf(Beta(k+a, n-k+b), p), 2ccdf(Beta(k+a, n-k+b), p))
end

# %%
function plot_pvalue_functions(; n = 20, k = 6, z = 2, a = 1, b = 1,
        f = trues(4), kwargs...)
    plot()
    f[1] && plot!(p -> pvalue_wilson_score(n, k, p), 0, 1;
        label="Wilson score", c=1)
    f[2] && plot!(p -> pvalue_bayesian(n, k, p; a, b), 0, 1;
        label="Bayesian", ls=:dashdot, c=2)
    f[3] && plot!(p -> pvalue_adjusted_wald(n, k, p; z), 0, 1;
        label=(z == 0 ? "" : "adjusted ") * "Wald", ls=:dash, c=3)
    f[4] && plot!(p -> pvalue_clopper_pearson(n, k, p), 0, 1;
        label="Clopper-Pearson", c=4)
    plot!(; xtick=0:0.1:1, ytick=0:0.1:1)
    title!("data: (n,k)=($n,$k)" * 
        (f[2] ? ", prior: Beta($a, $b)" : "") *
        (f[3] && z !== 0 ? ", adjustment for Wald: z=$z" : ""))
    plot!(; xguide="rate parameter θ", yguide="P-value")
    plot!(; kwargs...)
end

plot_pvalue_functions(; n = 20, k = 6, z = 2, a = 1, b = 1)

# %%
plot_pvalue_functions(; n = 20, k = 0, f = Bool[1,0,1,0], z = 0)

# %%
plot_pvalue_functions(; n = 20, k = 6, f = Bool[1,1,1,0], a=0.5, b=0.5)

# %%
using Roots
using StatsFuns
using DataFrames

# %%
function confint_clopper_pearson(n, k; α = 0.05)
    p_L = k > 0 ? quantile(Beta(k, n-k+1), α/2) : zero(α)
    p_U = k < n ? quantile(Beta(k+1, n-k), 1-α/2) : one(α)
    [p_L, p_U]
end

# %%
x ⪅ y = x < y || x ≈ y
_pdf_le(x, (dist, y)) =  pdf(dist, x) ⪅ y

function _search_boundary(f, x0, Δx, param)
    x = x0
    if f(x, param)
        while f(x - Δx, param) x -= Δx end
    else
        x += Δx
        while !f(x, param) x += Δx end
    end
    x
end

function pvalue_sterne(dist::DiscreteUnivariateDistribution, x)
    Px = pdf(dist, x)
    Px == 0 && return Px
    Px == 1 && return Px
    m = mode(dist)
    Px ≈ pdf(dist, m) && return one(Px)
    if x < m
        y = _search_boundary(_pdf_le, 2m - x, 1, (dist, Px))
        cdf(dist, x) + ccdf(dist, y-1)
    else # x > m
        y = _search_boundary(_pdf_le, 2m - x, -1, (dist, Px))
        cdf(dist, y) + ccdf(dist, x-1)
    end
end

function pvalue_sterne(n, k, p)
    pvalue_sterne(Binomial(n, p), k)
end

# 大きな n についてもうまく行くように
# Sterneの信頼区間の実装は難しい.
function confint_sterne(n, k; α = 0.05)
    a, b = confint_clopper_pearson(n, k; α = α/10)
    ps = find_zeros(a-√eps(), b+√eps()) do p
        logistic(0 < p ≤ 1 ? pvalue_sterne(n, k, p) : zero(p)) - logistic(α)
    end
    # 次の行は稀に区間にならない場合への対策
    [first(ps), last(ps)]
end

# %%
hcat(confint_sterne.(10, 0:10)...)'

# %%
function probabilities_of_type_I_error(; n = 100, z = 2, a = 1, b = 1, α = 0.05, L = 10^6)
    c_clopper_pearson = 0
    c_wilson_score = 0
    c_adjusted_wald = 0
    c_bayesian = 0
    c_sterne = 0
    for i in 1:L
        p = rand()
        k = rand(Binomial(n, p))
        c_clopper_pearson += pvalue_clopper_pearson(n, k, p) < α
        c_wilson_score += pvalue_wilson_score(n, k, p) < α
        c_adjusted_wald += pvalue_adjusted_wald(n, k, p; z) < α
        c_bayesian += pvalue_bayesian(n, k, p; a, b) < α
        c_sterne += pvalue_sterne(n, k, p) < α
    end
    DataFrame(
        method = [
            "Clopper-Pearson",
            "Sterne",
            "adjusted Wald",
            "Wilson score",
            "Bayesian",
        ],
        var"prob. of α-error" = [
            c_clopper_pearson,
            c_sterne,
            c_adjusted_wald,
            c_wilson_score,
            c_bayesian,
        ]/L,
        var"nominal α" = fill(α, 5),
        var"sample size" = fill(n, 5),
   )
end

# %%
@time probabilities_of_type_I_error(n = 5, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 10, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 20, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 30, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 50, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 100, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 200, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 300, α = 0.05)

# %%
@time probabilities_of_type_I_error(n = 1000, α = 0.05)

# %%
