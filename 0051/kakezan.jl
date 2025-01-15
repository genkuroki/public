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

# %% [markdown]
# https://blog.yuugakujyuku.com/2023/01/06/kakezantest/
#
# <img src="IMG_7744.jpeg" width=800>

# %%
using Distributions
using Printf
using Roots
using StatsFuns
using StatsPlots
default(fmt=:png)

# %%
safemul(x, y) = x == 0 ? x : isinf(x) ? oftype(x, Inf) : x*y
safediv(x, y) = x == 0 ? x : isinf(y) ? zero(y) : x/y

x ⪅ y = x < y || x ≈ y

mypdf(dist, x) = pdf(dist, x)
mypdf(dist::DiscreteUnivariateDistribution, x) = pdf(dist, round(Int, x))

distname(dist::Distribution) = replace(string(dist), r"{.*}" => "")
myskewness(dist) = skewness(dist)
mykurtosis(dist) = kurtosis(dist)
function standardized_moment(dist::ContinuousUnivariateDistribution, m)
    μ, σ = mean(dist), std(dist)
    quadgk(x -> (x - μ)^m * pdf(dist, x), extrema(dist)...)[1] / σ^m
end
myskewness(dist::MixtureModel{Univariate, Continuous}) =
    standardized_moment(dist, 3)
mykurtosis(dist::MixtureModel{Univariate, Continuous}) =
    standardized_moment(dist, 4) - 3

# %%
_riskratiohat(a, b, c, d) = safediv(a*(c+d), (a+b)*c)

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

function pvalue_rr_pearson_chisq(a, b, c, d; ρ=1.0)
    χ² = chisqstat_rr(a, b, c, d; ρ)
    ccdf(Chisq(1), χ²)
end

function confint_rr_pearson_chisq(a, b, c, d; α=0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0.0, Inf]
    f(logρ) = logit(pvalue_rr_pearson_chisq(a, b, c, d; ρ=exp(logρ))) - logit(α)
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
### score P-value for rate difference

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
    abs(Δ) == 1 && return Inf
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
    g(Δ) = chisqstat_rd_score(a, b, c, d; Δ, alg) - χ²_α
    Δ0 = riskdiffhat_score(a, b, c, d)
    L = if g(-1 + eps()) > 0
        find_zero(g, (-1 + eps(), Δ0), alg)
    else
        -1.0
    end
    U = if g(1 - eps()) > 0
        find_zero(g, (Δ0, 1 - eps()), alg)
    else
        1.0
    end
    [L, U]
end

# %%
_data = [
    75 64
    40 72
    75 72
    73 80
    59 72
    19 40
    66 72
    38 44
]
mean(_data; dims=1)

# %%
_data = 0.01 * _data .* [73 25]

# %%
_data = round.(Int, _data)

# %%
data = [_data[:,2] (25 .- _data[:,2]) _data[:,1] (73 .- _data[:,1])]

# %%
for (i, A) in enumerate(eachrow(data))
    print("Problem $i:  data = ", A)
    @printf(",  P-value of RR=1 = %.2f%%", 100pvalue_rr_pearson_chisq(A...))
    @printf(",  95%% CI of RR = [%.3f, %.3f]\n", confint_rr_pearson_chisq(A...)...)
end

# %%
PP = []
for (i, A) in enumerate(eachrow(data))
    a, b, c, d = A
    R1 = round(Int, 100a/(a+b))
    R2 = round(Int, 100c/(c+d))
    P = plot(ρ ->pvalue_rr_pearson_chisq(A...; ρ), 0.3, 6; label="", c=i)
    vline!([1]; label="", c=:black)
    plot!(xscale=:log10)
    xtick = Any[0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2, 2.4, 3, 4, 5, 6]
    xtick =(xtick, string.(xtick))
    plot!(; xtick, ytick=0:0.05:1)
    plot!(xguide="correct answer rate ratio RR", 
        yguide="P-value\ncompatibility of data and RR")
    title!("Problem $i: $c/$(c+d)≈$(R2)% vs. $a/$(a+b)≈$(R1)%")
    push!(PP, P)
end

plot((isodd(i) ? PP[i÷2+1] : PP[i÷2+4] for i in 1:8)...; size=(1200, 1600), layout=(4, 2))
plot!(leftmargin=10Plots.mm)
plot!(guidefontsize=14)

# %%
PP = []
for (i, A) in enumerate(eachrow(data))
    a, b, c, d = A
    R1 = round(Int, 100a/(a+b))
    R2 = round(Int, 100c/(c+d))
    ρs = 10.0 .^ range(log10(0.3), log10(6), 1000)
    P = plot(ρs, ρ ->-log2(pvalue_rr_pearson_chisq(A...; ρ)); label="", c=i)
    vline!([1]; label="", c=:black)
    plot!(xscale=:log10)
    xtick = Any[0.3, 0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2, 2.4, 3, 4, 5, 6]
    xtick =(xtick, string.(xtick))
    plot!(; xtick)
    plot!(ylim=(-0.3, 10.3), ytick=0:10)
    plot!(xguide="correct answer rate ratio RR", 
        yguide="S-value = −log\${}_2\$ P-value [bit]")
    title!("Problem $i: $c/$(c+d)≈$(R2)% vs. $a/$(a+b)≈$(R1)%")
    push!(PP, P)
end

plot((isodd(i) ? PP[i÷2+1] : PP[i÷2+4] for i in 1:8)...; size=(1200, 1600), layout=(4, 2))
plot!(leftmargin=12Plots.mm)
plot!(guidefontsize=14)

# %%
for (i, A) in enumerate(eachrow(data))
    print("Problem $i:  data = ", A)
    @printf(",  P-value of RD=0 = %.2f%%", 100pvalue_rd_score(A...))
    @printf(",  95%% CI of RD = [%.3f, %.3f]\n", confint_rd_score(A...)...)
end

# %%
PP = []
for (i, A) in enumerate(eachrow(data))
    a, b, c, d = A
    R1 = round(Int, 100a/(a+b))
    R2 = round(Int, 100c/(c+d))
    P = plot(Δ ->pvalue_rd_score(A...; Δ), -0.55, 0.65; label="", c=i)
    vline!([0]; label="", c=:black)
    plot!(; xtick=-1:0.1:1, ytick=0:0.05:1)
    plot!(xguide="correct answer rate difference RD", 
        yguide="P-value\ncompatibility of data and RD")
    title!("Problem $i: $c/$(c+d)≈$(R2)% vs. $a/$(a+b)≈$(R1)%")
    push!(PP, P)
end

plot((isodd(i) ? PP[i÷2+1] : PP[i÷2+4] for i in 1:8)...; size=(1200, 1600), layout=(4, 2))
plot!(leftmargin=10Plots.mm)
plot!(guidefontsize=14)

# %%
PP = []
for (i, A) in enumerate(eachrow(data))
    a, b, c, d = A
    R1 = round(Int, 100a/(a+b))
    R2 = round(Int, 100c/(c+d))
    P = plot(Δ -> -log2(pvalue_rd_score(A...; Δ)), -0.55, 0.65; label="", c=i)
    vline!([0]; label="", c=:black)
    plot!(; xtick=-1:0.1:1)
    plot!(; ylim=(-0.3, 10.3), ytick=0:10)
    plot!(xguide="correct answer rate difference RD", 
        yguide="S-value = −log\$_2\$ P-value [bit]")
    title!("Problem $i: $c/$(c+d)≈$(R2)% vs. $a/$(a+b)≈$(R1)%")
    push!(PP, P)
end

plot((isodd(i) ? PP[i÷2+1] : PP[i÷2+4] for i in 1:8)...; size=(1200, 1600), layout=(4, 2))
plot!(leftmargin=10Plots.mm)
plot!(guidefontsize=14)

# %%
