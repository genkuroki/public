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
using HypothesisTests
using Optim
using QuadGK
using Random
using Roots
using StatsFuns
using StatsPlots
default(fmt=:png, legendfontsize=10)

pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE=[raw"\usepackage{luatexja}"]
@eval Plots pgfx_sanitize_string(s::AbstractString) = s
gr()

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
        #ORhat = oddsratiohat(a, b, c, d)
        #ω_L, ω_U = ORhat/2, 2ORhat
        #[exp(find_zero(f, log(ω_L))), exp(find_zero(f, log(ω_U)))]
        logORhat = log(oddsratiohat(a, b, c, d))
        logL = find_zero(f, (logORhat-100, logORhat))
        logU = find_zero(f, (logORhat, logORhat+100))
        [exp(logL), exp(logU)]
    end
end

# %%
# Fisher's method (central) for OR

function pvalue_or_fisher_central(a, b, c, d; ω=1, correction=0.0)
    fnch = if ω == 1
        Hypergeometric(a+b, c+d, a+c)
    else
        FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)
    end
    min(1, 2(cdf(fnch, a) - correction*pdf(fnch, a)), 2(ccdf(fnch, a-1) - correction*pdf(fnch, a)))
end

function confint_or_fisher_central(a, b, c, d; α = 0.05, correction=0.0)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0, Inf]
    f(ω) = pvalue_or_fisher_central(a, b, c, d; ω, correction) - α
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

function pvalue_fisher_minlike(dist::DiscreteUnivariateDistribution, x; correction=0.0)
    Px = pdf(dist, x)
    Px == 0 && return Px
    Px == 1 && return Px
    m = mode(dist)
    if Px ≈ pdf(dist, m)
        (x ≠ m-1 && x ≠ m+1) ? one(Px) : 1 - pdf(dist, x)
    elseif x < m
        y = _search_boundary(_pdf_le, 2m - x, 1, (dist, Px))
            (cdf(dist, x) - correction*pdf(dist, x)) + (ccdf(dist, y-1) - correction*pdf(dist, y))
    else # x > m
        y = _search_boundary(_pdf_le, 2m - x, -1, (dist, Px))
        (cdf(dist, y) - correction*pdf(dist, y)) + (ccdf(dist, x-1) - correction*pdf(dist, x))
    end
end

function pvalue_or_fisher_minlike(a, b, c, d; ω=1, correction=0.0)
    fnch = if ω == 1
        Hypergeometric(a+b, c+d, a+c)
    else
        FisherNoncentralHypergeometric(a+b, c+d, a+c, ω)
    end
    pvalue_fisher_minlike(fnch, a; correction)
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
function pvals(randabcd; niters=10^5)
    pval_score = zeros(niters)
    pval_central_midp = zeros(niters)
    pval_minlike_midp = zeros(niters)
    Threads.@threads :static for i in 1:niters
        a, b, c, d = randabcd()
        pval_score[i] = pvalue_or_score(a, b, c, d)
        pval_central_midp[i] = pvalue_or_fisher_central(a, b, c, d; correction=0.5)
        pval_minlike_midp[i] = pvalue_or_fisher_minlike(a, b, c, d; correction=0.5)
    end
    (; pval_score, pval_central_midp, pval_minlike_midp)
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
    (; pval_score, pval_central_midp, pval_minlike_midp) = pvals(randabcd_2binomials; niters)
    
    print("probability of P-value ≤ 5%")
    println(p == q ? "  (α-error rate)" : "  (power)")
    println("  score:         ", round(100_ecdf(pval_score, 0.05); digits=1), "%")
    println("  minlike mid-p: ", round(100_ecdf(pval_minlike_midp, 0.05); digits=1), "%")
    println("  central mid-p: ", round(100_ecdf(pval_central_midp, 0.05); digits=1), "%")
    println()
    
    plot(α -> _ecdf(pval_score, α), 0, 0.1; label="score test")
    plot!(α -> _ecdf(pval_minlike_midp, α), 0, 0.1; label="minlike mid-p", ls=:dash)
    plot!(α -> _ecdf(pval_central_midp, α), 0, 0.1; label="central mid-p", ls=:dot)
    p == q && plot!(identity; label="", ls=:dot, c=:black, alpha=0.5)
    plot!(; legend)
    plot!(; xtick=0:0.01:1, ytick, xrotation=90)
    plot!(; xguide="α", yguide)
    title!("model: $(distname(bin1))×$(distname(bin2))", titlefontsize=12)
    plot!(size=(400, 400))
end

# %%
m, n = 10, 10
for q in 0.1:0.1:0.5
    plot_pvals(; bin1=Binomial(m, q), bin2=Binomial(n, q), niters=10^6) |> display
end

# %%
m, n = 8, 12
for q in 0.1:0.1:0.5
    plot_pvals(; bin1=Binomial(m, q), bin2=Binomial(n, q), niters=10^6) |> display
end

# %%
m, n = 5, 10
for q in 0.1:0.1:0.5
    plot_pvals(; bin1=Binomial(m, q), bin2=Binomial(n, q), niters=10^6) |> display
end

# %%
m, n = 16, 20
for q in 0.05:0.05:0.5
    plot_pvals(; bin1=Binomial(m, q), bin2=Binomial(n, q), niters=10^6) |> display
end

# %%
a, b, c, d = 1, 10, 5, 5

ωmin, ωmax = 0.001, 5
ωs = exp.(range(log(ωmin), log(ωmax), 500))
plot(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.0); label="score")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.2); label="score (0.2-corrected)", ls=:dot)
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.5); label="score (0.5-orrected)", ls=:dot)
#plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.0); label="Fisher (minlik)")
plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.5); label="Fisher (minlike, mid-p)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.0); label="Fisher (central)")
plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.5); label="Fisher (central, mid-p)", ls=:dot)
xtick = Any[0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
xtick = xtick[ωmin .≤ xtick .≤ ωmax]
xtick = (xtick, string.(xtick))
plot!(; xscale=:log, xtick, ytick=0:0.05:1, xrotation=90)
plot!(xguide="OR", yguide="P-value")
plot!(legend=:topleft)
title!("data: $([a b; c d])")

# %%
a, b, c, d = 1, 100, 4, 50

ωmin, ωmax = 0.001, 5
ωs = exp.(range(log(ωmin), log(ωmax), 500))
plot(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.0); label="score")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.2); label="score (0.2-corrected)")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.5); label="score (0.5-orrected)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.0); label="Fisher (minlik)")
plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.5); label="Fisher (minlike, mid-p)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.0); label="Fisher (central)")
plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.5); label="Fisher (central, mid-p)", ls=:dot)
xtick = Any[0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
xtick = xtick[ωmin .≤ xtick .≤ ωmax]
xtick = (xtick, string.(xtick))
plot!(; xscale=:log, xtick, ytick=0:0.05:1, xrotation=90)
plot!(xguide="OR", yguide="P-value")
plot!(legend=:topleft)
title!("data: $([a b; c d])")

# %%
a, b, c, d = 1, 24, 5, 15

ωmin, ωmax = 0.001, 5
ωs = exp.(range(log(ωmin), log(ωmax), 500))
plot(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.0); label="score")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.2); label="score (0.2-corrected)")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.5); label="score (0.5-orrected)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.0); label="Fisher (minlik)")
plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.5); label="Fisher (minlike, mid-p)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.0); label="Fisher (central)")
plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.5); label="Fisher (central, mid-p)", ls=:dot)
xtick = Any[0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
xtick = xtick[ωmin .≤ xtick .≤ ωmax]
xtick = (xtick, string.(xtick))
plot!(; xscale=:log, xtick, ytick=0:0.05:1, xrotation=90)
plot!(xguide="OR", yguide="P-value")
plot!(legend=:topleft)
title!("data: $([a b; c d])")

# %%
a, b, c, d = 4, 3107-4, 67, 10883-67

ωmin, ωmax = 0.02, 2
ωs = exp.(range(log(ωmin), log(ωmax), 500))
plot(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.0); label="score")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.2); label="score (0.2-corrected)")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.5); label="score (0.5-orrected)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.0); label="Fisher (minlik)")
plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.5); label="Fisher (minlike, mid-p)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.0); label="Fisher (central)")
plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.5); label="Fisher (central, mid-p)", ls=:dot)
xtick = Any[0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
xtick = xtick[ωmin .≤ xtick .≤ ωmax]
xtick = (xtick, string.(xtick))
plot!(; xscale=:log, xtick, ytick=0:0.05:1, xrotation=90)
plot!(legend=:topleft)
plot!(xguide="OR", yguide="P-value")
plot!(legend=:topleft)
title!("data: $([a b; c d])")

# %%
a, b, c, d = 4, 3107-4, 67, 10883-67

ωmin, ωmax = 0.2, 1
ωs = exp.(range(log(ωmin), log(ωmax), 500))
plot(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.0); label="score")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.2); label="score (0.2-corrected)")
#plot!(ωs, ω -> pvalue_or_score(a, b, c, d; ω, correction=0.5); label="score (0.5-orrected)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.0); label="Fisher (minlik)")
plot!(ωs, ω -> pvalue_or_fisher_minlike(a, b, c, d; ω, correction=0.5); label="Fisher (minlike, mid-p)", ls=:dash)
#plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.0); label="Fisher (central)")
plot!(ωs, ω -> pvalue_or_fisher_central(a, b, c, d; ω, correction=0.5); label="Fisher (central, mid-p)", ls=:dot)
xtick = Any[0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
xtick = xtick[ωmin .≤ xtick .≤ ωmax]
xtick = (xtick, string.(xtick))
plot!(; xscale=:log, xtick, ytick=0:0.05:1, xrotation=90)
plot!(legend=:topleft)
plot!(xguide="OR", yguide="P-value")
plot!(legend=:topleft)
title!("data: $([a b; c d])")

# %%
