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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Roots
using StatsFuns
using StatsPlots
default(fmt=:png)

x ⪅ y = x < y || x ≈ y
_ecdf(A, x) = count(≤(x), A) / length(A)

pvalue_central(dist, x) = min(1, 2cdf(dist, x), 2ccdf(dist, x-1))

function pvalue_minlike(dist, x)
    m = mode(dist)
    x == m && return 1.0
    px = pdf(dist, x)
    if x < m
        y = m
        while !(pdf(dist, y) ⪅ px) y += 1 end
        cdf(dist, x) + ccdf(dist, y-1)
    else
        y = m
        while !(pdf(dist, y) ⪅ px) y -= 1 end
        cdf(dist, y) + ccdf(dist, x-1)
    end
end

function pvalue_normal_approx(dist, x; correction=0.0)
    absz = max(0, abs(x - mean(dist)) - correction) / std(dist)
    2ccdf(Normal(), absz)
end

function pvalue_gamma_approx_eti(dist, x; correction=0.0)
    μ, σ² = mean(dist), var(dist)
    κ, θ = μ^2/σ², σ²/μ
    gamma = Gamma(κ, θ)
    m = median(gamma)
    min(1, 2cdf(gamma, min(m, x + correction)), 2ccdf(gamma, max(m, x - correction)))
end

function pvalue_hdi(dist, x;
        xmin = minimum(dist) == -Inf ? -1e8 : minimum(dist),
        xmax = maximum(dist) ==  Inf ?  1e8 : maximum(dist),
        correction = 0.0,
    )
    m = mode(dist)
    px = pdf(dist, x)
    px ≈ pdf(dist, m) && return 1.0
    f(y) = pdf(dist, y) - px
    if x < m
        y = find_zero(f, (m, 1.0))
        cdf(dist, min(m, x + correction)) + ccdf(dist, max(m, y - correction))
    else
        y = find_zero(f, (0.0, m))
        cdf(dist, min(m, y + correction)) + ccdf(dist, max(m, x - correction))
    end
end

function pvalue_gamma_approx_hdi(dist, x; correction=0.0)
    μ, σ² = mean(dist), var(dist)
    κ, θ = μ^2/σ², σ²/μ
    gamma = Gamma(κ, θ)
    pvalue_hdi(gamma, x; correction)
end

loglikrat(dist1, dist0, x) = 2(logpdf(dist1, x) - logpdf(dist0, x))
pvalue_loglikrat(dist1, dist0, x; df = 1) = ccdf(Chisq(df), loglikrat(dist1, dist0, x))

function pvalue_score(k, n, p)
    p̂ = k/n
    χ² = (p̂ - p)^2 / (p*(1-p)/n)
    ccdf(Chisq(1), χ²)
end

function pvalue_logmaxlikrat(k, n, p)
    p̂ = k/n
    χ² = 2(xlogy(k, p̂/p) + xlogy(n-k, (1-p̂)/(1-p)))
    ccdf(Chisq(1), χ²)
end

function pvalue_bayes_eti(k, n, p; prior=Beta(1, 1))
    κ, λ = params(prior)
    beta = Beta(κ+k, λ+n-k)
    min(1, 2cdf(beta, p), 2ccdf(beta, p))
end

function pvalue_bayes_hdi(k, n, p; prior=Beta(1, 1))
    κ, λ = params(prior)
    beta = Beta(κ+k, λ+n-k)
    pvalue_hdi(beta, p)
end

function sim_bin(n, p; L=10^5)
    bin = Binomial(n, p)
    pval_score = zeros(L)
    pval_lmlr = zeros(L)
    pval_central = zeros(L)
    pval_minlike = zeros(L)
    Threads.@threads for i in 1:L
        k = rand(bin)
        pval_score[i] = pvalue_score(k, n, p)
        pval_lmlr[i] = pvalue_logmaxlikrat(k, n, p)
        pval_central[i] = pvalue_central(bin, k)
        pval_minlike[i] = pvalue_minlike(bin, k)
    end
    pval_score, pval_lmlr, pval_central, pval_minlike
end

function sim_negbin(k, p; L=10^5)
    negbin = NegativeBinomial(k, p) + k
    pval_score = zeros(L)
    pval_lmlr = zeros(L)
    pval_central = zeros(L)
    pval_minlike = zeros(L)
    Threads.@threads for i in 1:L
        n = rand(negbin)
        pval_score[i] = pvalue_score(k, n, p)
        pval_lmlr[i] = pvalue_logmaxlikrat(k, n, p)
        pval_central[i] = pvalue_central(negbin, n)
        pval_minlike[i] = pvalue_minlike(negbin, n)
    end
    pval_score, pval_lmlr, pval_central, pval_minlike
end

function pdf_invnegbin(x; k=7, p=0.5) # x = k/n
    0 < x ≤ 1 || return 0.0
    negbin = NegativeBinomial(k, p) + k
    n = round(Int, k/x)
    w = k/max(k, n-1/2) - k/(n+1/2)
    pdf(negbin, n) / w
end

function plot_invnegbin(; k=7, p=0.5, xmin=0.0, xmax=1.0)
    negbin = NegativeBinomial(k, p) + k
    @show μ = sum(k/n * pdf(negbin, n) for n in k:10^3)
    @show σ = √sum((k/n - μ)^2 * pdf(negbin, n) for n in k:10^3)
    @show k/μ, mean(negbin)
    plot(x -> pdf_invnegbin(x; k, p), xmin, xmax; label="inv.neg.bin.")
    plot!(Normal(μ, σ), xmin, xmax; label="normal approx.")
end

function print_result(k, n, p; prior=Beta(1, 1))
    @show k n p
    bin = bin0 = Binomial(n, p)
    bin1 = Binomial(n, k/n)
    negbin = negbin0 = NegativeBinomial(k, p) + k
    negbin1 = NegativeBinomial(k, k/n) + k
    println()
    #@show mean(bin) std(bin)
    #@show (k - mean(bin)) / std(bin)
    #@show loglikrat(bin1, bin0, k)
    #@show pvalue_loglikrat(bin1, bin0, k)
    @show pvalue_normal_approx(bin, k)
    @show pvalue_central(bin, k)
    @show pvalue_normal_approx(bin, k; correction=0.5)
    @show pvalue_minlike(bin, k)
    println()
    #@show mean(negbin) std(negbin)
    #@show (n - mean(negbin)) / std(negbin)
    #@show loglikrat(negbin1, negbin0, n)
    #@show pvalue_loglikrat(negbin1, negbin0, n)
    @show pvalue_central(negbin, n)
    @show pvalue_gamma_approx_eti(negbin, n; correction=0.5)
    #@show pvalue_gamma_approx_eti(negbin, n)
    @show pvalue_minlike(negbin, n)
    #@show pvalue_normal_approx(negbin, n; correction=0.5)
    #@show pvalue_normal_approx(negbin, n)
    @show pvalue_gamma_approx_hdi(negbin, n; correction=0.5)
    #@show pvalue_gamma_approx_hdi(negbin, n)
    println()
    #@show pvalue_logmaxlikrat(k, n, p)
    @show pvalue_score(k, n, p)
    #@show pvalue_bayes_eti(k, n, p; prior)
    @show pvalue_bayes_hdi(k, n, p; prior)
end

function plot_negbin(k, p)
    negbin = NegativeBinomial(k, p) + k
    #@show negbin
    @show μ, σ² = mean(negbin), var(negbin)
    @show κ, θ = μ^2/σ², σ²/μ
    σ = √σ²
    P1 = bar(negbin; alpha=0.3, label="NegativeBinomial($k, $p) + $k")
    plot!(Normal(μ, σ); label="normal approx.", lw=2)
    plot!(Gamma(κ, θ); label="gamma approx.", ls=:dash, lw=2)
    plot!(legend=:outertop, xlim=(μ-5σ, μ+5σ))
end

function print_and_plot_result(k, n, p; prior=Beta(1, 1))
    print_result(k, n, p; prior)
    println()
    plot_negbin(k, p)
end

# %%
plot_invnegbin(k=7, p=0.5)

# %%
print_and_plot_result(12, 24, 0.3)

# %%
print_and_plot_result(33, 80, 0.3)

# %%
print_and_plot_result(87, 240, 0.3)

# %%
print_and_plot_result(7, 24, 1/2)

# %%
print_and_plot_result(31, 80, 1/2)

# %%
print_and_plot_result(105, 240, 1/2)

# %%
print_and_plot_result(469, 1000, 1/2)

# %%
pval_score, pval_lmlr, pval_central, pval_minlike = sim_bin(24, 0.5)

plot(α -> _ecdf(pval_score, α), 0, 0.16; label="score")
#plot!(α -> _ecdf(pval_lmlr, α); label="log max.lik.rat.", ls=:dash, c=2)
#plot!(α -> _ecdf(pval_minlike, α); label="minlike", ls=:dashdotdot, c=3)
plot!(α -> _ecdf(pval_central, α); label="central", ls=:dashdot, c=4)
plot!(identity; label="", c=:black, ls=:dot, lw=0.7)
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=90)
plot!(size=(400, 400))

# %%
pval_score, pval_lmlr, pval_central, pval_minlike = sim_bin(50, 0.25)

plot(α -> _ecdf(pval_score, α), 0, 0.1; label="score")
#plot!(α -> _ecdf(pval_lmlr, α); label="log max.lik.rat.", ls=:dash, c=2)
#plot!(α -> _ecdf(pval_minlike, α); label="minlike", ls=:dashdotdot, c=3)
plot!(α -> _ecdf(pval_central, α); label="central", ls=:dashdot, c=4)
plot!(identity; label="", c=:black, ls=:dot, lw=0.7)
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=90)
plot!(size=(400, 400))

# %%
pval_score, pval_lmlr, pval_central, pval_minlike = sim_bin(50, 0.4)

plot(α -> _ecdf(pval_score, α), 0, 0.1; label="score")
plot!(α -> _ecdf(pval_lmlr, α); label="log max.lik.rat.", ls=:dash, c=2)
plot!(α -> _ecdf(pval_minlike, α); label="minlike", ls=:dashdotdot, c=3)
plot!(α -> _ecdf(pval_central, α); label="central", ls=:dashdot, c=4)
plot!(identity; label="", c=:black, ls=:dot, lw=0.7)
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=90)
plot!(size=(400, 400))

# %%
pval_score, pval_lmlr, pval_central, pval_minlike = sim_negbin(7, 0.5)

plot(α -> _ecdf(pval_score, α), 0, 0.16; label="score")
#plot!(α -> _ecdf(pval_lmlr, α); label="log max.lik.rat.", ls=:dash, c=2)
#plot!(α -> _ecdf(pval_minlike, α); label="minlike", ls=:dashdotdot, c=3)
plot!(α -> _ecdf(pval_central, α); label="central", ls=:dashdot, c=4)
plot!(identity; label="", c=:black, ls=:dot, lw=0.7)
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=90)
plot!(size=(400, 400))

# %%
pval_score, pval_lmlr, pval_central, pval_minlike = sim_negbin(7, 0.5)

plot(α -> _ecdf(pval_score, α), 0, 0.16; label="score")
#plot!(α -> _ecdf(pval_lmlr, α); label="log max.lik.rat.", ls=:dash, c=2)
plot!(α -> _ecdf(pval_minlike, α); label="minlike", ls=:dashdotdot, c=3)
#plot!(α -> _ecdf(pval_central, α); label="central", ls=:dashdot, c=4)o
plot!(identity; label="", c=:black, ls=:dot, lw=0.7)
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=90)
plot!(size=(400, 400))

# %%
pval_score, pval_lmlr, pval_central, pval_minlike = sim_negbin(10, 0.6)

plot(α -> _ecdf(pval_score, α), 0, 0.1; label="score")
plot!(α -> _ecdf(pval_lmlr, α); label="log max.lik.rat.", ls=:dash, c=2)
plot!(α -> _ecdf(pval_minlike, α); label="minlike", ls=:dashdotdot, c=3)
plot!(α -> _ecdf(pval_central, α); label="central", ls=:dashdot, c=4)
plot!(identity; label="", c=:black, ls=:dot, lw=0.7)
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=90)
plot!(size=(400, 400))

# %%
function plot_negbin_chisq(; k=7, p=0.5, xmin=0.0, xmax=8.0, L=10^5)
    negbin = NegativeBinomial(k, p) + k
    scorechisq = [(k/n - p)^2/(p*(1-p)/n) for n in (rand(negbin) for _ in 1:L)]
    plot(x -> _ecdf(scorechisq, x), xmin, xmax; label="ecdf of score χ² statistics")
    plot!(x -> cdf(Chisq(1), x); label="cdf of Chisq(1)")
end

plot_negbin_chisq(k=7, p=0.5)

# %%
plot_negbin_chisq(k=20, p=0.5)

# %%
function plot_negbin_z(; k=7, p=0.5, xmin=-4.5, xmax=4.5, L=10^5)
    negbin = NegativeBinomial(k, p) + k
    scorez = [(k/n - p)/√(p*(1-p)/n) for n in (rand(negbin) for _ in 1:L)]
    plot(x -> _ecdf(scorez, x), xmin, xmax; label="ecdf of score z statistics")
    plot!(x -> cdf(Normal(), x), xmin, xmax; label="cdf of Normal()")
end

plot_negbin_z(k=7, p=0.5)

# %%
plot_negbin_z(k=40, p=0.2)

# %%
