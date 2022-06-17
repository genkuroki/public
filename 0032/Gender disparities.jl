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

# %% [markdown]
# https://twitter.com/drk0311/status/1537076913788751872
#
# https://www.resuscitationjournal.com/article/S0300-9572(20)30105-2/fulltext
#
# ![FVTKldcVUAAWbjX.jpg](attachment:b4865a57-14fc-4065-b009-9061d2ea6e17.jpg)

# %%
ENV["LINES"], ENV["COLUMNS"] = 100, 300
using Distributions
using StatsPlots
default(fmt=:png, size=(400, 250), titlefontsize=10)
using KernelDensity
using DataFrames
using Roots
using StatsFuns

safemul(x, y) = x == 0 ? x : isinf(x) ? typeof(x)(Inf) : x*y
safediv(x, y) = x == 0 ? x : isinf(y) ? zero(y) : x/y

# %%
function pvalue_wilson(n, k, p)
    p̂ = k/n
    SE = √(p*(1-p)/n)
    2ccdf(Normal(), abs(p̂ - p)/SE)
end

function confint_wilson(n, k; α = 0.05)
    p̂ = k/n
    z = quantile(Normal(), 1-α/2)
    a, b, c = 1+z^2/n, p̂+z^2/(2n), p̂^2
    # ap² - 2bp + c = 0 を解く.
    sqrtD = √(b^2 - a*c)
    p_L = (b - sqrtD)/a
    p_U = (b + sqrtD)/a
    [p_L, p_U]
end

# %%
oddsratiohat(a, b, c, d) = safediv(a*d, b*c)
stderr_logoddsratiohat(a, b, c, d) = √(1/a + 1/b + 1/c + 1/d)

function pvalue_or_wald(a, b, c, d; ω=1)
    logORhat = log(oddsratiohat(a, b, c, d))
    SEhat_logORhat = stderr_logoddsratiohat(a, b, c, d)
    2ccdf(Normal(0, 1), safediv(abs(logORhat - log(ω)), SEhat_logORhat))
end

function confint_or_wald(a, b, c, d; α=0.05)
    z = quantile(Normal(), 1-α/2)
    ORhat = oddsratiohat(a, b, c, d)
    SEhat_logORhat = stderr_logoddsratiohat(a, b, c, d)
    [safemul(exp(-z*SEhat_logORhat), ORhat), safemul(exp(z*SEhat_logORhat), ORhat)]
end

# %%
riskratiohat(a, b, c, d) = safediv(a*(c+d), (a+b)*c)
stderr_logriskratiohat(a, b, c, d) = √(1/a - 1/(a+b) + 1/c - 1/(c+d))

function pvalue_rr_wald(a, b, c, d; ρ=1)
    logRRhat = log(riskratiohat(a, b, c, d))
    SEhat_logRRhat = stderr_logriskratiohat(a, b, c, d)
    2ccdf(Normal(0, 1), safediv(abs(logRRhat - log(ρ)), SEhat_logRRhat))
end

function confint_rr_wald(a, b, c, d; α=0.05)
    z = quantile(Normal(), 1-α/2)
    RRhat = riskratiohat(a, b, c, d)
    SEhat_logRRhat = stderr_logriskratiohat(a, b, c, d)
    [safemul(exp(-z*SEhat_logRRhat), RRhat), safemul(exp(z*SEhat_logRRhat), RRhat)]
end

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

function pvalue_or_pearson(a, b, c, d; ω=1, correction=0.0)
    χ² = chisqstat_or(a, b, c, d; ω, correction)
    ccdf(Chisq(1), χ²)
end

function confint_or_pearson(a, b, c, d; α=0.05, correction=0.0)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0.0, Inf]
    f(ω) = logit(pvalue_or_pearson(a, b, c, d; ω, correction)) - logit(α)
    if a == 0 || d == 0
        [0.0, find_zero(f, 1.0)]
    elseif b == 0 || c == 0
        [find_zero(f, 1.0), Inf]
    else
        ω_L, ω_U = confint_or_wald(a, b, c, d; α)
        [find_zero(f, ω_L), find_zero(f, ω_U)]
    end
end

# %%
function Delta(a, b, c, d; ρ=1)
    m, n = a+b, c+d
    A, B, C = ρ-1, n-a+ρ*(m-c), a*n-ρ*m*c
    isinf(ρ) ? typeof(ω)(-c) : safediv(2C, B + √(B^2 - 4A*C))
end

function _chisqstat_rr(a, b, c, d, Δ)
    m, n = a+b, c+d
    safemul(Δ^2, safediv(b, m*(a-Δ)) + safediv(d, n*(c+Δ)))
end

function chisqstat_rr(a, b, c, d; ρ=1)
    Δ = Delta(a, b, c, d; ρ)
    _chisqstat_rr(a, b, c, d, Δ)
end

function pvalue_rr_pearson(a, b, c, d; ρ=1)
    χ² = chisqstat_rr(a, b, c, d; ρ)
    ccdf(Chisq(1), χ²)
end

function confint_rr_pearson(a, b, c, d; α=0.05)
    (a+b==0 || c+d==0 || a+c==0 || b+d==0) && return [0.0, Inf]
    f(ρ) = logit(pvalue_rr_pearson(a, b, c, d; ρ)) - logit(α)
    if a == 0 || d == 0
        [0.0, find_zero(f, 1.0)]
    elseif b == 0 || c == 0
        [find_zero(f, 1.0), Inf]
    else
        ρ_L, ρ_U = confint_rr_wald(a, b, c, d; α)
        [find_zero(f, ρ_L), find_zero(f, ρ_U)]
    end
end

# %%
aed = [
      0   5   0   3
      1  14   1  13
      9 101   2  53
     14 172   2  76
     41 305   6 125
     88 529  16 142
    144 762  23 150
    124 857  22 206
     63 460  32 206
      7 108   9  71
]

age = [
     0  4
     5 14
    15 24
    25 34
    35 44
    45 54
    55 64
    65 74
    75 84
    85 120
]

name = [
    "min of age"
    "max of age"
    "Male AED"
    "Male data size"
    "Female AED"
    "Female data size"
]

df = DataFrame([age aed], name)

# %%
Male_CI = vcat((confint_wilson(n, k)' for (n, k) in zip(df."Male data size", df."Male AED"))...)
Female_CI = vcat((confint_wilson(n, k)' for (n, k) in zip(df."Female data size", df."Female AED"))...)
df."Male CI L" = Male_CI[:,1]
df."Male CI U" = Male_CI[:,2]
df."Female CI L" = Female_CI[:,1]
df."Female CI U" = Female_CI[:,2]
df

# %%
xticklabel = ["""$a-$(b > 100 ? "" : b)""" for (a,b) in zip(df."min of age", df."max of age")]
xtick = (1:10, xticklabel)
a = df."Male AED"
m = df."Male data size"
M = @. a/m
MCI_L = df."Male CI L"
MCI_U = df."Male CI U"
c = df."Female AED"
n = df."Female data size"
F = @. c/n
FCI_L = df."Female CI L"
FCI_U = df."Female CI U"

plot()
plot!(M; label="Male", marker=:o, ribbon=(M - MCI_L, MCI_U - M))
plot!(F; label="Female", marker=:o, ribbon=(F - FCI_L, FCI_U - F))
plot!(; xtick)

# %%
ORhat = @. oddsratiohat(a, m-a, c, n-c)

# %%
α = 0.05
ORCI = vcat((confint_or_pearson(a, m-a, c, n-c; α)' for (m, a, n, c)
        in zip(df."Male data size", df."Male AED", df."Female data size", df."Female AED"))...)
ORCI_L = ORCI[:,1]
ORCI_U = ORCI[:,2]

plot(ORhat; label="odds ratio odds(M)/odds(F)", ribbon=(ORhat-ORCI_L, ORCI_U-ORhat))
hline!([1]; label="")
plot!(; xtick, ytick=0:20)
title!("$(100(1-α))% confidence interval of odds ratio")

# %%
α = 0.1
ORCI = vcat((confint_or_pearson(a, m-a, c, n-c; α)' for (m, a, n, c)
        in zip(df."Male data size", df."Male AED", df."Female data size", df."Female AED"))...)
ORCI_L = ORCI[:,1]
ORCI_U = ORCI[:,2]

plot(ORhat; label="odds ratio odds(M)/odds(F)", ribbon=(ORhat-ORCI_L, ORCI_U-ORhat))
hline!([1]; label="")
plot!(; xtick, ytick=0:20)
title!("$(100(1-α))% confidence interval of odds ratio")

# %%
RRhat = @. riskratiohat(a, m-a, c, n-c)

# %%
α = 0.05
RRCI = vcat((confint_rr_pearson(a, m-a, c, n-c; α)' for (m, a, n, c)
        in zip(df."Male data size", df."Male AED", df."Female data size", df."Female AED"))...)
RRCI_L = RRCI[:,1]
RRCI_U = RRCI[:,2]

plot(RRhat; label="percentage ratio M/F", ribbon=(RRhat-RRCI_L, RRCI_U-RRhat))
hline!([1]; label="")
plot!(; xtick, ytick=0:20)
title!("$(100(1-α))% confidence interval of percentage ratio")

# %%
α = 0.1
RRCI = vcat((confint_rr_pearson(a, m-a, c, n-c; α)' for (m, a, n, c)
        in zip(df."Male data size", df."Male AED", df."Female data size", df."Female AED"))...)
RRCI_L = RRCI[:,1]
RRCI_U = RRCI[:,2]

plot(RRhat; label="percentage ratio M/F", ribbon=(RRhat-RRCI_L, RRCI_U-RRhat))
hline!([1]; label="")
plot!(; xtick, ytick=0:20)
title!("$(100(1-α))% confidence interval of percentage ratio")

# %%
function plot_pvalues(; idx = 3:8, xlim = (1/15, 15), df=df)
    @show agemin = df."min of age"[idx[begin]]
    @show agemax = df."max of age"[idx[end]]
    @show a = sum(df."Male AED"[idx])
    @show m = sum(df."Male data size"[idx])
    @show c = sum(df."Female AED"[idx])
    @show n = sum(df."Female data size"[idx])
    println()

    @show ORhat = oddsratiohat(a, m-a, c, n-c)
    @show RRhat = riskratiohat(a, m-a, c, n-c)
    @show pvalue_or_pearson(a, m-a, c, n-c)
    @show pvalue_rr_pearson(a, m-a, c, n-c)
    @show confint_or_pearson(a, m-a, c, n-c)
    @show confint_rr_pearson(a, m-a, c, n-c)

    xtick_ = vcat((10.0^k*[1, 2, 5] for k in -1:2)...)
    xtick_ = xtick_[@.(xlim[begin] ≤ xtick_ ≤ xlim[end])]
    @show xtick_
    xtick = (xtick_, string.(xtick_))
    
    P = plot(ω -> pvalue_or_pearson(a, m-a, c, n-c; ω), xlim...; label="")
    vline!([ORhat]; label="$(round(ORhat; digits=2))", ls=:dash, c=1)
    vline!([1]; label="", c=2)
    plot!(; xscale=:log10, xtick, ytick=0:0.1:1)
    plot!(; xguide="hypothetical odds ratio M/F", yguide="P-value")
    plot!(; bottommargin=4Plots.mm, leftmargin=4Plots.mm)
    
    Q = plot(ρ -> pvalue_rr_pearson(a, m-a, c, n-c; ρ), xlim...; label="")
    vline!([RRhat]; label="$(round(RRhat; digits=2))", ls=:dash, c=1)
    vline!([1]; label="", c=2)
    plot!(; xscale=:log10, xtick, ytick=0:0.1:1)
    plot!(; xguide="hypothetical percentage ratio M/F", yguide="P-value")
    plot!(; bottommargin=4Plots.mm, leftmargin=4Plots.mm)
    
    Title = plot(; framestyle=:none, bottommargin=-50Plots.px)
    title!("""Age: $agemin - $(agemax > 100 ? "infinity" : agemax),   Data: Male $a/$m,  Female $c/$n""")

    plot(Title, P, Q; size=(800, 250), layout=@layout[a{0.1h}; [b c]])
    plot!(; titlefontsize=9, tickfontsize=6, guidefontsize=8)
end

# %%
P_all = plot_pvalues(idx = 1:10)

# %%
P15_74 = plot_pvalues(idx = 3:8)

# %%
P1 = plot_pvalues(idx = 1:2)

# %%
P2 = plot_pvalues(idx = 3:4)

# %%
P3 = plot_pvalues(idx = 5:6)

# %%
P4 = plot_pvalues(idx = 7:8)

# %%
P5 = plot_pvalues(idx = 9:10)

# %%
plot(P1, P2, P3, P4, P5; size=(800, 1250), layout=(5, 1))

# %%
aed2 = aed[begin:2:end, :] + aed[begin+1:2:end, :]
age2 = vcat(([age[i,1] age[i+1,2]] for i in 1:2:10)...)
df2 = DataFrame([age2 aed2], name)

Male2_CI = vcat((confint_wilson(n, k)' for (n, k) in zip(df2."Male data size", df2."Male AED"))...)
Female2_CI = vcat((confint_wilson(n, k)' for (n, k) in zip(df2."Female data size", df2."Female AED"))...)
df2."Male CI L" = Male2_CI[:,1]
df2."Male CI U" = Male2_CI[:,2]
df2."Female CI L" = Female2_CI[:,1]
df2."Female CI U" = Female2_CI[:,2]
df2

# %%
xticklabel2 = ["""$a-$(b > 100 ? "" : b)""" for (a,b) in zip(df2."min of age", df2."max of age")]
xtick2 = (1:5, xticklabel2)
a2 = df2."Male AED"
m2 = df2."Male data size"
M2 = @. a2/m2
MCI2_L = df2."Male CI L"
MCI2_U = df2."Male CI U"
c2 = df2."Female AED"
n2 = df2."Female data size"
F2 = @. c2/n2
FCI2_L = df2."Female CI L"
FCI2_U = df2."Female CI U"

plot()
plot!(M2; label="Male", marker=:o, ribbon=(M2 - MCI2_L, MCI2_U - M2))
plot!(F2; label="Female", marker=:o, ribbon=(F2 - FCI2_L, FCI2_U - F2))
plot!(; xtick=xtick2)

# %%
ORhat2 = @. oddsratiohat(a2, m2-a2, c2, n2-c2)

# %%
α = 0.05
ORCI2 = vcat((confint_or_pearson(a2, m2-a2, c2, n2-c2; α)' for (m2, a2, n2, c2)
        in zip(df2."Male data size", df2."Male AED", df2."Female data size", df2."Female AED"))...)
ORCI2_L = ORCI2[:,1]
ORCI2_U = ORCI2[:,2]

plot(ORhat2; label="odds ratio odds(M)/odds(F)", ribbon=(ORhat2-ORCI2_L, ORCI2_U-ORhat2))
hline!([1]; label="")
plot!(; xtick=xtick2, ytick=0:20)
title!("$(100(1-α))% confidence interval of odds ratio")

# %%
α = 0.10
ORCI2 = vcat((confint_or_pearson(a2, m2-a2, c2, n2-c2; α)' for (m2, a2, n2, c2)
        in zip(df2."Male data size", df2."Male AED", df2."Female data size", df2."Female AED"))...)
ORCI2_L = ORCI2[:,1]
ORCI2_U = ORCI2[:,2]

plot(ORhat2; label="odds ratio odds(M)/odds(F)", ribbon=(ORhat2-ORCI2_L, ORCI2_U-ORhat2))
hline!([1]; label="")
plot!(; xtick=xtick2, ytick=0:20)
title!("$(100(1-α))% confidence interval of odds ratio")

# %%
RRhat2 = @. riskratiohat(a2, m2-a2, c2, n2-c2)

# %%
α = 0.05
RRCI2 = vcat((confint_rr_pearson(a2, m2-a2, c2, n2-c2; α)' for (m2, a2, n2, c2)
        in zip(df2."Male data size", df2."Male AED", df2."Female data size", df2."Female AED"))...)
RRCI2_L = RRCI2[:,1]
RRCI2_U = RRCI2[:,2]

plot(RRhat2; label="percentage ratio M/F", ribbon=(RRhat2-RRCI2_L, RRCI2_U-RRhat2))
hline!([1]; label="")
plot!(; xtick=xtick2, ytick=0:20)
title!("$(100(1-α))% confidence interval of percentage ratio")

# %%
α = 0.10
RRCI2 = vcat((confint_rr_pearson(a2, m2-a2, c2, n2-c2; α)' for (m2, a2, n2, c2)
        in zip(df2."Male data size", df2."Male AED", df2."Female data size", df2."Female AED"))...)
RRCI2_L = RRCI2[:,1]
RRCI2_U = RRCI2[:,2]

plot(RRhat2; label="percentage ratio M/F", ribbon=(RRhat2-RRCI2_L, RRCI2_U-RRhat2))
hline!([1]; label="")
plot!(; xtick=xtick2, ytick=0:20)
title!("$(100(1-α))% confidence interval of percentage ratio")

# %%
