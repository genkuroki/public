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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions
using Roots
using StatsBase
using StatsPlots
default(fmt=:png, titlefontsize=8, tickfontsize=6, size=(400, 250),
    plot_titlefontsize=10)
safemul(x, y) = x == 0 ? x : isinf(x) ? typeof(x)(Inf) : x*y
safediv(x, y) = x == 0 ? x : isinf(y) ? zero(y) : x/y

# %%
# 警告: 以下の実装の精度は低い. 改良の仕方が
# 
# Kenneth J. Rothman, Sander Greenland, and Timothy L. Lash
# Modern Epistemology, Third Edition, 2008, 888 pages
#
# の
#
# Chapter 14. Instroduction to Categorical Statistics
# Section. Two Study Groups: Large-Sample Methods, pp.299-300
#
# に書いてある. そこでは, 次の文献が引用されている:
#
# Guangyong Zou and Allan Donner
# A simple alternative confidence interval for the difference between two proportions
# Controlled Clinical Trials, Volume 25, Issue 1, February 2004, Pages 3-12
# https://doi.org/10.1016/j.cct.2003.08.010
#
# Zou-Donnerの信頼区間に対応するP値函数の実装については
#
# https://github.com/genkuroki/public/blob/main/0033/probability%20of%20alpha%20error%20of%20Zou-Donner.ipynb
#
# を参照せよ.

riskdiffhat(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function stderr_riskdiffhat(a, b, c, d; u=0)
    m, n = a+b, c+d
    p̂, q̂ = safediv(a, m), safediv(c, n)
    √(safediv(p̂*(1-p̂), m-u) + safediv(q̂*(1-q̂), n-u))
end

function pvalue_rd_wald(a, b, c, d; Δ=0, u=0)
    RDhat = riskdiffhat(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat(a, b, c, d; u)
    2ccdf(Normal(0, 1), safediv(abs(RDhat - Δ), SEhat_riskdiffhat))
end

function confint_rd_wald(a, b, c, d; α=0.05, u=0)
    z = quantile(Normal(), 1-α/2)
    RDhat = riskdiffhat(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat(a, b, c, d; u)
    [RDhat - z*SEhat_riskdiffhat, RDhat + z*SEhat_riskdiffhat]
end

# %%
@show confint_rd_wald(0, 1, 1, 1) confint_rd_wald(1e-4, 1, 1, 1)
@show confint_rd_wald(1, 0, 1, 1) confint_rd_wald(1, 1e-4, 1, 1)
@show confint_rd_wald(1, 1, 0, 1) confint_rd_wald(1, 1, 1e-4, 1)
@show confint_rd_wald(1, 1, 1, 0) confint_rd_wald(1, 1, 1, 1e-4)
@show confint_rd_wald(0, 0, 1, 1) confint_rd_wald(1e-4, 1e-8, 1, 1)
@show confint_rd_wald(1, 1, 0, 0) confint_rd_wald(1, 1, 1e-4, 1e-8)
@show confint_rd_wald(0, 1, 0, 1) confint_rd_wald(1e-4, 1, 1e-8, 1)
@show confint_rd_wald(1, 0, 1, 0) confint_rd_wald(1, 1e-4, 1, 1e-8)
@show confint_rd_wald(0, 1, 1, 0) confint_rd_wald(1e-4, 1, 1, 1e-8)
@show confint_rd_wald(1, 0, 0, 1) confint_rd_wald(1, 1e-4, 1e-8, 1);

# %%
# risk difference Zou-Donner

riskdiffhat_zou_donner(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function stderr_riskdiffhat_zou_donner(a, b, c, d; u=1)
    m, n = a+b, c+d
    p̂, q̂ = safediv(a, m), safediv(c, n)
    √(safediv(p̂*(1-p̂), m-u) + safediv(q̂*(1-q̂), n-u))
end

function pvalue_rd_zou_donner(a, b, c, d; Δ=0, u=1)
    ((a==0 && d==0) || (b==0 && c==0)) && return 1.0
    RDhat = riskdiffhat_zou_donner(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat_zou_donner(a, b, c, d; u)
    Z = safediv((1 - RDhat^2)*abs(atanh(RDhat) - atanh(Δ)), SEhat_riskdiffhat)
    2ccdf(Normal(), abs(Z))
end

function confint_rd_zou_donner(a, b, c, d; α=0.05, u=1)
    z = quantile(Normal(), 1-α/2)
    RDhat = riskdiffhat_zou_donner(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat_zou_donner(a, b, c, d; u)
    m = atanh(RDhat)
    d = safediv(z*SEhat_riskdiffhat, 1 - RDhat^2)
    [tanh(m-d), tanh(m+d)]
end

# %%
@show confint_rd_zou_donner(0, 1, 1, 1) confint_rd_zou_donner(1e-4, 1, 1, 1)
@show confint_rd_zou_donner(1, 0, 1, 1) confint_rd_zou_donner(1, 1e-4, 1, 1)
@show confint_rd_zou_donner(1, 1, 0, 1) confint_rd_zou_donner(1, 1, 1e-4, 1)
@show confint_rd_zou_donner(1, 1, 1, 0) confint_rd_zou_donner(1, 1, 1, 1e-4)
@show confint_rd_zou_donner(0, 0, 1, 1) confint_rd_zou_donner(1e-4, 1e-8, 1, 1)
@show confint_rd_zou_donner(1, 1, 0, 0) confint_rd_zou_donner(1, 1, 1e-4, 1e-8)
@show confint_rd_zou_donner(0, 1, 0, 1) confint_rd_zou_donner(1e-4, 1, 1e-8, 1)
@show confint_rd_zou_donner(1, 0, 1, 0) confint_rd_zou_donner(1, 1e-4, 1, 1e-8)
@show confint_rd_zou_donner(0, 1, 1, 0) confint_rd_wald(1e-4, 1, 1, 1e-8)
@show confint_rd_zou_donner(1, 0, 0, 1) confint_rd_wald(1, 1e-4, 1e-8, 1);

# %%
function ci_rd(pvaluefunc; α=0.05)
    f(x) = pvaluefunc(x) - α
    find_zeros(f, -1, 1)
end

a, b, c, d = 426-255, 255, 433-277, 277

@show confint_rd_zou_donner(a, b, c, d; α=0.05)
f(Δ) = pvalue_rd_zou_donner(a, b, c, d; Δ)
@show ci_rd(f; α=0.05);

# %%
function pvalue_rd_zou_donner_prime(a, b, c, d; Δ=0)
    ((a==0 && d==0) || (b==0 && c==0)) && return 1.0
    RDhat = riskdiffhat(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat(a, b, c, d)
    Z = safediv((1 - Δ^2)*abs(atanh(RDhat) - atanh(Δ)), SEhat_riskdiffhat)
    2ccdf(Normal(), abs(Z))
end

# %%
@show a, b, c, d = 1, 9, 19, 1
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_zou_donner_prime(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)

plot(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ), -1, 1; label="")
#plot!(Δ -> pvalue_rd_zou_donner_prime(a, b, c, d; Δ), -1, 1; label="", ls=:dash)
plot!(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -1, 1; label="", ls=:dashdot)

# %%
@show a, b, c, d = 1, 9, 19, 1
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_zou_donner_prime(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)

plot(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ), -1, 1; label="")
plot!(Δ -> pvalue_rd_zou_donner_prime(a, b, c, d; Δ), -1, 1; label="", ls=:dash)
#plot!(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -1, 1; label="", ls=:dashdot)

# %%
function delta(a, b, c, d; ω=1)
    A, B, C = 1-ω, a+d+ω*(b+c), a*d-ω*b*c
    isinf(ω) ? typeof(ω)(-min(b, c)) : safediv(2C, B + √(B^2 - 4A*C))
end

# correction = 0.5 は連続性補正を与える.
function _chisqstat_or(a, b, c, d, δ; correction=0.0)
    ã, b̃, c̃, d̃ = a-δ, b+δ, c+δ, d-δ
    safemul(max(0, abs(δ)-correction)^2, 1/ã + 1/b̃ + 1/c̃ + 1/d̃)
end

function chisqstat_or(a, b, c, d; ω=1, correction=0.0)
    δ = delta(a, b, c, d; ω)
    _chisqstat_or(a, b, c, d, δ; correction)
end

function pvalue_or_pearson_chisq(a, b, c, d; ω=1, correction=0.0)
    χ² = chisqstat_or(a, b, c, d; ω, correction)
    ccdf(Chisq(1), χ²)
end

function confint_or_pearson_chisq(a, b, c, d; α=0.05, correction=0.0)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0, Inf]
    f(logω) = logit(pvalue_or_pearson_chisq(a, b, c, d; ω=exp(logω), correction)) - logit(α)
    ps = if a == 0 || d == 0
        [0, exp(find_zero(f, 0.0))]
    elseif b == 0 || c == 0
        [exp(find_zero(f, 0.0)), Inf]
    else
        ORhat = oddsratiohat(a, b, c, d)
        ω_L, ω_U = ORhat/2, 2ORhat
        [exp(find_zero(f, log(ω_L))), exp(find_zero(f, log(ω_U)))]
    end
end

# %%
"""Bayes版P値函数達を作る函数"""
function make_pvalue_rd_bayes(a, b, c, d; M=10^6, conjprior=(0.5, 0.5))
    α, β = conjprior
    p = rand(Beta(α+a, β+b), M)
    q = rand(Beta(α+c, β+d), M)
    Δ = @. p - q
    ecdf_RD = ecdf(Δ)
    pvalue_rd_bayes(Δ) = min(1, 2ecdf_RD(Δ), 2(1-ecdf_RD(Δ)))
    pvalue_rd_bayes
end

# %%
a, b, c, d = 10, 3, 4, 9
pvalue_rd_bayes = make_pvalue_rd_bayes(a, b, c, d)
plot(legend=:topleft)
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ), -1, 1; label="ZD")
plot!(Δ -> pvalue_rd_wald(a, b, c, d; Δ); label="Wald", ls=:dash)
plot!(Δ -> pvalue_rd_bayes(Δ); label="Bayes", ls=:dashdotdot)

# %%
oddsratiohat(a, b, c, d) = safediv(a*d, b*c)

function sim_alphaerrors(m, n, p, q=p; L=10^5)
    Δ = p - q
    ω = safediv(p*(1-q), (1-p)*q)
    bin1, bin2 = Binomial(m, p), Binomial(n, q)
    pval_wald = similar(zeros(), L)
    pval_zou_donner = similar(zeros(), L)
    pval_chisq = similar(zeros(), L)
    Threads.@threads for i in 1:L
        a, c = rand(bin1), rand(bin2)
        b, d = m-a, n-c
        pval_wald[i] = pvalue_rd_wald(a, b, c, d; Δ)
        pval_zou_donner[i] = pvalue_rd_zou_donner(a, b, c, d; Δ)
        pval_chisq[i] = pvalue_or_pearson_chisq(a, b, c, d; ω)
    end
    ecdf_wald = ecdf(pval_wald)
    ecdf_zou_donner = ecdf(pval_zou_donner)
    ecdf_chisq = ecdf(pval_chisq)
    F_wald(x) = ecdf_wald(x)
    F_zou_donner(x) = ecdf_zou_donner(x)
    F_chisq(x) = ecdf_chisq(x)
    (; F_wald, F_zou_donner, F_chisq)
end

function plot_alphaerrors(m, n, p, q=p; L=10^5, kwargs...)
    (; F_wald, F_zou_donner, F_chisq) = sim_alphaerrors(m, n, p, q)
    plot(legend=:bottomright)
    plot!(F_zou_donner, 0, 0.1; label="ZD")
    plot!(F_wald, 0, 0.1; label="Wald", ls=:dash)
    plot!(F_chisq, 0, 0.1; label="chisq", ls=:dashdot)
    plot!(identity, 0, 0.1; label="", c=:red, ls=:dot, alpha=0.7)
    plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
    title!("Bin(m=$m, p=$p)×Bin(n=$n, q=$q)")
    plot!(size=(400, 400))
    plot!(; kwargs...)
end

function plot_alphaerrors3x3(m, n; L=10^5, kwargs...)
    for p in 0.1:0.1:0.5
        PP = []
        for q in 0.1:0.1:0.9
            P = plot_alphaerrors(m, n, p, q; L, kwargs...)
            push!(PP, P)
        end
        plot(PP...; size=(800, 840), layout=(3, 3))
        plot!(plot_title="m = $m, n = $n, p = $p") |> display
    end
end

# %%
plot_alphaerrors3x3(10, 20)

# %%
plot_alphaerrors3x3(20, 40)

# %%
plot_alphaerrors3x3(30, 60)

# %%
plot_alphaerrors3x3(40, 80)

# %%
plot_alphaerrors3x3(10, 10)

# %%
plot_alphaerrors3x3(20, 20)

# %%
plot_alphaerrors3x3(30, 30)

# %%
plot_alphaerrors3x3(40, 40)

# %%
