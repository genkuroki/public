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
using QuadGK
using Random
using Roots
using StatsPlots
default(fmt=:png, legendfontsize=11)

distname(dist) = replace(string(dist),
    r"{[^}]*}"=>"",
    r".="=>"",
    "Binomial"=>"Bin"
)

_ecdf(A, x) = count(≤(x), A) / length(A)

safediv(x, y) = x == 0 ? zero(x/y) : x/y

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
    plot!(size=(500, 500))
end

function plot_powers(; )
    
end

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.25), bin2=Binomial(8, 0.25))

# %% tags=[]
plot_pvals(; bin1=Binomial(6, 0.25), bin2=Binomial(8, 0.89))

# %% tags=[]
plot_pvals(; bin1=Binomial(12, 0.25), bin2=Binomial(16, 0.25))

# %%
plot_pvals(; bin1=Binomial(12, 0.25), bin2=Binomial(16, 0.755))

# %% tags=[]
plot_pvals(; bin1=Binomial(24, 0.25), bin2=Binomial(32, 0.25))

# %% tags=[]
plot_pvals(; bin1=Binomial(24, 0.25), bin2=Binomial(32, 0.615))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.25), bin2=Binomial(160, 0.25))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.25), bin2=Binomial(160, 0.41))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.025), bin2=Binomial(160, 0.025))

# %% tags=[]
plot_pvals(; bin1=Binomial(120, 0.025), bin2=Binomial(160, 0.105))

# %%
