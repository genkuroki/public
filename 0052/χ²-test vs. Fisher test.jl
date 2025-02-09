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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %%
using Distributions
using HypothesisTests
using Optim
using QuadGK
using Random
using Roots
using StatsFuns
using StatsPlots
default(fmt=:png, legendfontsize=10)

distname(dist) = replace(string(dist), r"{[^}]*}"=>"", r".="=>"")

function distname(dist::Binomial)
    n, p = params(dist)
    p = round(p; sigdigits=3)
    "Bin($n, $p)"
end

safemul(x, y) = x == 0 ? zero(x/y) : isinf(x) ? oftype(x, Inf) : x*y
safediv(x, y) = x == 0 ? zero(x/y) : isinf(y) ? zero(y) : x/y
x ⪅ y = x < y || x ≈ y
_ecdf(A, x) = count(≤(x), A) / length(A)

function pvalue_chisq(a, b, c, d; yates=false)
    N = a+b+c+d
    chisqstat = safediv(N*max(0, abs(a*d - b*c) - (N/2)yates)^2, (a+b)*(c+d)*(a+c)*(b+d))
    ccdf(Chisq(1), chisqstat)
end

function pvals(randabcd; niters=10^5)
    pval_chisq = zeros(niters)
    pval_yates = zeros(niters)
    pval_minlike = zeros(niters)
    pval_central = zeros(niters)
    Threads.@threads :static for i in 1:niters
        a, b, c, d = randabcd()
        pval_chisq[i] = pvalue_chisq(a, b, c, d)
        pval_yates[i] = pvalue_chisq(a, b, c, d; yates=true)
        fet = FisherExactTest(a, b, c, d)
        pval_minlike[i] = pvalue(fet; method=:minlike)
        pval_central[i] = pvalue(fet)
    end
    (; pval_chisq, pval_yates, pval_minlike, pval_central)
end

function make_randabcd(bin1::Binomial, bin2::Binomial)
    function randabcd_2binomials()
        a, c = rand(bin1), rand(bin2)
        b, d = ntrials(bin1)-a, ntrials(bin2)-c
        a, b, c, d
    end
    randabcd_2binomials
end

function plot_pvals(; bin1=Binomial(12, 0.25), bin2=Binomial(16, 0.25), niters=10^6)
    m, p = params(bin1)
    n, q = params(bin2)
    ytick = p == q ? (0:0.01:1) : (0:0.05:1)
    legend = p == q ? true : (:bottomright)
    yguide = "probability of P-value ≤ α"
    yguide *= p == q ? "  (α-error rate)" : "  (power)"
    
    println("true expectation value: ")
    Base.print_array(stdout, round.([m*p m*(1-p); n*q n*(1-q)]; sigdigits=3))
    println("\n")
    
    randabcd_2binomials = make_randabcd(bin1, bin2)
    (; pval_chisq, pval_yates, pval_minlike, pval_central) = pvals(randabcd_2binomials; niters)
    
    print("probability of P-value ≤ 5%")
    println(p == q ? "  (α-error rate)" : "  (power)")
    println("  χ²-test:               ", round(100_ecdf(pval_chisq, 0.05); digits=1), "%")
    println("  Fisher test (minlike): ", round(100_ecdf(pval_minlike, 0.05); digits=1), "%")
    println("  χ²-test (Yates):       ", round(100_ecdf(pval_yates, 0.05); digits=1), "%")
    println("  Fisher test (central): ", round(100_ecdf(pval_central, 0.05); digits=1), "%")
    println()
    
    plot(α -> _ecdf(pval_chisq, α), 0, 0.1; label="χ²-test")
    plot!(α -> _ecdf(pval_minlike, α), 0, 0.1; label="Fisher test (minlike)", ls=:dash)
    plot!(α -> _ecdf(pval_yates, α), 0, 0.1; label="χ²-test (Yates)", ls=:dot)
    plot!(α -> _ecdf(pval_central, α), 0, 0.1; label="Fisher test (central)", ls=:dashdot)
    p == q && plot!(identity; label="", ls=:dot, c=:black, alpha=0.5)
    plot!(; legend)
    plot!(; xtick=0:0.01:1, ytick, xrotation=90)
    plot!(; xguide="α", yguide)
    title!("model: $(distname(bin1))×$(distname(bin2))", titlefontsize=12)
    plot!(size=(400, 400))
end

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.25), bin2=Binomial(8, 0.25))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.25), bin2=Binomial(8, 0.8))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.25), bin2=Binomial(8, 0.89))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.1), bin2=Binomial(8, 0.1))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.2), bin2=Binomial(8, 0.2))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.3), bin2=Binomial(8, 0.3))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.4), bin2=Binomial(8, 0.4))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.5), bin2=Binomial(8, 0.5))

# %% tags=[]
plot_pvals(; bin1=Binomial(10, 0.1), bin2=Binomial(10, 0.1))

# %% tags=[]
plot_pvals(; bin1=Binomial(10, 0.2), bin2=Binomial(10, 0.2))

# %% tags=[]
plot_pvals(; bin1=Binomial(10, 0.3), bin2=Binomial(10, 0.3))

# %% tags=[]
plot_pvals(; bin1=Binomial(10, 0.4), bin2=Binomial(10, 0.4))

# %% tags=[]
plot_pvals(; bin1=Binomial(10, 0.5), bin2=Binomial(10, 0.5))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.25), bin2=Binomial(16, 0.25))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.25), bin2=Binomial(16, 0.65))

# %%
plot_pvals(; bin1=Binomial(12, 0.25), bin2=Binomial(16, 0.755))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.1), bin2=Binomial(16, 0.1))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.2), bin2=Binomial(16, 0.2))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.3), bin2=Binomial(16, 0.3))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.4), bin2=Binomial(16, 0.4))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.5), bin2=Binomial(16, 0.5))

# %% tags=[]
plot_pvals(; bin1=Binomial(24, 0.25), bin2=Binomial(32, 0.25))

# %% tags=[]
plot_pvals(; bin1=Binomial(24, 0.25), bin2=Binomial(32, 0.55))

# %% tags=[]
plot_pvals(; bin1=Binomial(24, 0.25), bin2=Binomial(32, 0.615))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.25), bin2=Binomial(160, 0.25))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.25), bin2=Binomial(160, 0.38))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.25), bin2=Binomial(160, 0.41))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.025), bin2=Binomial(160, 0.025))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.025), bin2=Binomial(160, 0.09))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.025), bin2=Binomial(160, 0.105))

# %% tags=[]
plot_pvals(; bin1=Binomial(1200, 0.0025), bin2=Binomial(1600, 0.0025))

# %% tags=[]
plot_pvals(; bin1=Binomial(1200, 0.0025), bin2=Binomial(1600, 0.0095))

# %% tags=[]
plot_pvals(; bin1=Binomial(1200, 0.0025), bin2=Binomial(1600, 0.0112))

# %% tags=[]
plot_pvals(; bin1=Binomial(2400, 0.001), bin2=Binomial(3200, 0.001))

# %% tags=[]
plot_pvals(; bin1=Binomial(2400, 0.001), bin2=Binomial(3200, 0.004))

# %% tags=[]
plot_pvals(; bin1=Binomial(2400, 0.001), bin2=Binomial(3200, 0.0051))

# %%
oddsratiohat(a, b, c, d) = safediv(a*d, b*c)

# score method for OR

function delta(a, b, c, d; ω=1)
    A, B, C = 1-ω, a+d+ω*(b+c), a*d-ω*b*c
    isinf(ω) ? oftype(ω, -min(b, c)) : safediv(2C, B + √(B^2 - 4A*C))
end

# correction = 0.5 は連続性補正を与える.
function _chisqstat_or(a, b, c, d, δ; correction=0.0)
    ã, b̃, c̃, d̃ = a-δ, b+δ, c+δ, d-δ
    safemul(max(0, abs(δ)-correction)^2, 1/ã + 1/b̃ + 1/c̃ + 1/d̃)
end

function chisqstat_or(a, b, c, d; ω=1, correction=0.0)
    δ = delta(a, b, c, d; ω)
    _chisqstat_or(a, b, c, d, δ; correction)
end

function pvalue_or_score(a, b, c, d; ω=1, correction=0.0)
    χ² = chisqstat_or(a, b, c, d; ω, correction)
    ccdf(Chisq(1), χ²)
end

function confint_or_score(a, b, c, d; α=0.05, correction=0.0)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0, Inf]
    f(logω) = pvalue_or_score(a, b, c, d; ω=exp(logω), correction) - α
    ps = if a == 0 || d == 0
        [0, exp(find_zero(f, 0.0))]
    elseif b == 0 || c == 0
        [exp(find_zero(f, 0.0)), Inf]
    else
        ORhat = oddsratiohat(a, b, c, d)
        ω_L, ω_U = ORhat/10, 10ORhat
        [exp(find_zero(f, log(ω_L))), exp(find_zero(f, log(ω_U)))]
    end
end

# %%
# score method for RR

_riskratiohat(a, b, c, d) = safediv(a*(c+d), (a+b)*c)

# ((a-Δ)/(a-Δ+b))/((c+Δ)/(c+Δ+d)) = ρ if Δ = Delta(a, b, c, d; ρ)
function Delta(a, b, c, d; ρ=1.0)
    m, n = a+b, c+d
    A, B, C = ρ-1, n-a+ρ*(m-c), a*n-ρ*m*c
    Δ = isinf(ρ) ? oftype(ρ, -c) : ρ==0 ? oftype(ρ, a) : safediv(2C, B + √(B^2 - 4A*C))
end

function _chisqstat_rr(a, b, c, d, Δ)
    m, n = a+b, c+d
    safemul(Δ^2, safediv(b, m*(a-Δ)) + safediv(d, n*(c+Δ)))
end

function chisqstat_rr(a, b, c, d; ρ=1.0)
    Δ = Delta(a, b, c, d; ρ)
    _chisqstat_rr(a, b, c, d, Δ)
end

function pvalue_rr_score(a, b, c, d; ρ=1.0)
    χ² = chisqstat_rr(a, b, c, d; ρ)
    ccdf(Chisq(1), χ²)
end

function confint_rr_score(a, b, c, d; α=0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0.0, Inf]
    f(logρ) = logit(pvalue_rr_score(a, b, c, d; ρ=exp(logρ))) - logit(α)
    L = if f(-Inf) > 0
        -Inf
    else
        logRRhat = log(_riskratiohat(a, b, c, d))
        x0 = logRRhat == -Inf ? -10.0 : logRRhat == Inf ? 10.0 : logRRhat - 1
        find_zero(f, x0)
    end
    U = if f(Inf) > 0
        Inf
    else
        logRRhat = log(_riskratiohat(a, b, c, d))
        x0 = logRRhat == -Inf ? -10.0 : logRRhat == Inf ? 10.0 : logRRhat + 1
        find_zero(f, x0)
    end
    [exp(L), exp(U)]
end

# %%
### score method for RD

riskdiffhat_score(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function loglik_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safemul(a, log(p)) + safemul(b, log(1-p)) + safemul(c, log(q)) + safemul(d, log(1-q))
end

function scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(a, p) - safediv(b, 1-p) + safediv(c, q) - safediv(d, 1-q)
end

function d_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    -safediv(a, p^2) - safediv(b, (1-p)^2) - safediv(c, q^2) - safediv(d, (1-q)^2)
end

function scorestat_Δ_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(a, p) - safediv(b, 1-p)
end

function estimate_q_given_Δ_rd(a, b, c, d, Δ=0.0; alg=Bisection())
    qmin, qmax = max(0.0, -Δ), min(1.0, 1.0-Δ)
    a+c==0 && return qmin
    b+d==0 && return qmax
    f(q) = scorestat_q_rd(a, b, c, d, q, Δ)
    S_qmin = f(qmin + eps())
    S_qmax = f(qmax - eps())
    S_qmin ≥ 0 && S_qmax ≥ 0 && return S_qmin < S_qmax ? qmin : qmax
    S_qmin ≤ 0 && S_qmax ≤ 0 && return S_qmin < S_qmax ? qmax : qmin
    find_zero(f, (qmin + eps(), qmax - eps()), alg)
end

function varinv_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(p*(1-p), a+b) + safediv(q*(1-q), c+d)
end

function chisqstat_rd_score(a, b, c, d; Δ=0.0, alg=Bisection())
    Δ = clamp(Δ, -1 + eps(), 1 - eps())
    q̃ = estimate_q_given_Δ_rd(a, b, c, d, Δ; alg)
    S = scorestat_Δ_rd(a, b, c, d, q̃, Δ)
    Vinv = varinv_scorestat_q_rd(a, b, c, d, q̃, Δ)
    safemul(S^2, Vinv)
end

function pvalue_rd_score(a, b, c, d; Δ=0.0, alg=Bisection())
    χ² = chisqstat_rd_score(a, b, c, d; Δ, alg)
    ccdf(Chisq(1), χ²)
end

function confint_rd_score(a, b, c, d; α=0.05, alg=Bisection())
    χ²_α = cquantile(Chisq(1), α)
    RDhat = riskdiffhat_score(a, b, c, d)
    g(Δ) = chisqstat_rd_score(a, b, c, d; Δ, alg) - χ²_α
    L = if g(-1 + eps()) > 0
        find_zero(g, (-1 + eps(), RDhat), alg)
    else
        -1.0
    end
    U = if g(1 - eps()) > 0
        find_zero(g, (RDhat, 1 - eps()), alg)
    else
        1.0
    end
    [L, U]
end

# %%
# Fisher's method (minlike) for OR

function oddsratiohat_fisher(a, b, c, d)
    f(x) = -logpdf(FisherNoncentralHypergeometric(a+b, c+d, a+c, x[1]), a)
    ω₀ = oddsratiohat(a, b, c, d)
    o = optimize(f, [ω₀])
    o.minimizer[1]
end

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

function pvalue_fisher_minlike(dist::DiscreteUnivariateDistribution, x)
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

function pvalue_or_fisher_minlike(a, b, c, d; ω=1)
    fnch = if ω == 1
        Hypergeometric(a+b, c+d, a+c)
    else
        FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)
    end
    pvalue_fisher_minlike(fnch, a)
end

function find_pos(f, x)
    while f(x) ≤ 0
        x *= 2
    end
    x
end

function confint_or_fisher_minlike(a, b, c, d; α = 0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0, Inf]
    f(logω) = logit(pvalue_or_fisher_minlike(a, b, c, d; ω=exp(logω))) - logit(α)
    if a == 0 || d == 0
        [0.0, exp(find_zero(f, (find_pos(f, -1.0), 0.0)))]
    elseif b == 0 || c == 0
        [exp(find_zero(f, (0.0, find_pos(f, 1.0)))), Inf]
    else
        ω_L, ω_U = confint_or_score(a, b, c, d; α)
        ω_L, ω_U = ω_L/(ω_U/ω_L), ω_U*(ω_U/ω_L)
        ps = exp.(find_zeros(f, log(ω_L), log(ω_U)))
        # 次の行は稀に区間にならない場合への対策
        isempty(ps) ? [0, Inf] : [first(ps), last(ps)]
    end
end

# %%
# Fisher's method (central) for OR

function pvalue_or_fisher_central(a, b, c, d; ω=1)
    fnch = if ω == 1
        Hypergeometric(a+b, c+d, a+c)
    else
        FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)
    end
    min(1, 2cdf(fnch, a), 2ccdf(fnch, a-1))
end

function confint_or_fisher_central(a, b, c, d; α = 0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0, Inf]
    f(ω) = logit(pvalue_or_fisher_central(a, b, c, d; ω)) - logit(α)
    if a == 0 || d == 0
        [0.0, find_zero(f, 1.0)]
    elseif b == 0 || c == 0
        [find_zero(f, 1.0), Inf]
    else
        ω_L, ω_U = confint_or_score(a, b, c, d; α)
        ω_L, ω_U = ω_L/(ω_U/ω_L), ω_U*(ω_U/ω_L)
        find_zeros(f, ω_L, ω_U)
    end
end

# %%
r(x) = round(x; sigdigits=3)

# %% [markdown]
# https://drmagician.exblog.jp/22086293/ より
#
# <img src="IMG_7998.png">

# %%
a, b, c, d = 1, 24, 5, 15
@show a, b, c, d
println()

@show pvalue_or_score(a, b, c, d; correction=0.0) |> r
@show pvalue_or_score(a, b, c, d; correction=0.5) |> r
@show pvalue_or_fisher_minlike(a, b, c, d) |> r
@show pvalue_or_fisher_central(a, b, c, d) |> r
println()

@show oddsratiohat(a, b, c, d)
@show oddsratiohat_fisher(a, b, c, d)
println()

@show confint_or_score(a, b, c, d; correction=0.0) .|> r
@show confint_or_score(a, b, c, d; correction=0.5) .|> r
@show confint_or_fisher_minlike(a, b, c, d) .|> r
@show confint_or_fisher_central(a, b, c, d) .|> r
println()

@show _riskratiohat(a, b, c, d)
@show confint_rr_score(a, b, c, d) .|> r
println()

@show riskdiffhat_score(a, b, c, d)
@show confint_rd_score(a, b, c, d) .|> r
;

# %% [markdown]
# https://drmagician.exblog.jp/22086293/ より
#
# <img src="IMG_7998.png">

# %%
@show a, b, c, d
@show oddsratiohat(a, b, c, d)
@show oddsratiohat_fisher(a, b, c, d)

ωmin, ωmax = 0.0008, 5
plot()
plot!(ω -> pvalue_or_score(a, b, c, d; ω, correction=0.0), ωmin, ωmax; label="score")
plot!(ω -> pvalue_or_fisher_minlike(a, b, c, d; ω), ωmin, ωmax; label="Fisher (minlike)", ls=:dash)
plot!(ω -> pvalue_or_score(a, b, c, d; ω, correction=0.5), ωmin, ωmax; label="score (Yates)", ls=:dot)
plot!(ω -> pvalue_or_fisher_central(a, b, c, d; ω), ωmin, ωmax; label="Fisher (central)", ls=:dashdot)
plot!(legend=:topleft, legendfontsize=12)
xtick = Any[0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]
xtick = (xtick, string.(xtick))
plot!(; xscale=:log, xtick, ytick=0:0.05:1)
plot!(xguide="hypothetical odds ratio", yguide="P-value", guidefontsize=14)
title!("data: $([a b; c d])")

# %% [markdown]
# https://drmagician.exblog.jp/22086293/ より
#
# <img src="IMG_7998.png">

# %%
@show a, b, c, d
@show _riskratiohat(a, b, c, d)

ρmin, ρmax = 0.005, 5
plot(ρ -> pvalue_rr_score(a, b, c, d; ρ), ρmin, ρmax; label="score")
plot!(legend=:topleft, legendfontsize=12)
xtick = Any[0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]
xtick = (xtick, string.(xtick))
plot!(; xscale=:log, xtick, ytick=0:0.05:1)
plot!(xguide="hypothetical risk ratio", yguide="P-value", guidefontsize=14)
title!("data: $([a b; c d])")

# %% [markdown]
# https://drmagician.exblog.jp/22086293/ より
#
# <img src="IMG_7998.png">

# %%
@show a, b, c, d
@show riskdiffhat_score(a, b, c, d)

plot(Δ -> pvalue_rd_score(a, b, c, d; Δ), -0.6, 0.2; label="score")
plot!(legend=:topleft, legendfontsize=12)
plot!(xtick=-1:0.1:1, ytick=0:0.05:1)
plot!(xguide="hypothetical risk difference", yguide="P-value", guidefontsize=14)
title!("data: $([a b; c d])")

# %% [markdown]
# https://drmagician.exblog.jp/22086293/ より
#
# <img src="IMG_7998.png"> <img src="IMG_7999.png">

# %%
m, n = 25, 20
q = 6/45
plot_pvals(; bin1=Binomial(m, q), bin2=Binomial(n, q))

# %% [markdown]
# https://drmagician.exblog.jp/22086293/ より
#
# <img src="IMG_7998.png"> <img src="IMG_7999.png">

# %%
m, n = 25, 20
p, q = 1/25, 5/20
@show p-q |> r
@show p/q |> r
@show (p/(1-p))/(q/(1-q)) |> r
println()
plot_pvals(; bin1=Binomial(m, p), bin2=Binomial(n, q))

# %%
