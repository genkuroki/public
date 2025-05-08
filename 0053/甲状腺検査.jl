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
#     display_name: Julia current stable release
#     language: julia
#     name: julia
# ---

# %%
using Distributions
using Roots
using StatsFuns
using StatsPlots
default(fmt=:png)

safemul(x, y) = x == 0 ? x : isinf(x) ? oftype(x, Inf) : x*y
safediv(x, y) = x == 0 ? x : isinf(y) ? zero(y) : x/y

oddsratiohat(a, b, c, d) = safediv(a*d, b*c)

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

# %% [markdown]
# <img src="IMG_8592.png">

# %%
A = [
    14 41981-14
     1 18327-1
]
a, b, c, d = A'

@show A
@show oddsratiohat(a, b, c, d)
@show confint_or_pearson_chisq(a, b, c, d)
@show pvalue_or_pearson_chisq(a, b, c, d)
;

# %%
A = [
    14 41981-14
     2 18327-2
]
a, b, c, d = A'

@show A
@show oddsratiohat(a, b, c, d)
@show confint_or_pearson_chisq(a, b, c, d)
@show pvalue_or_pearson_chisq(a, b, c, d)
;

# %%
A = [
    11 17969-11
     1 18327-1
]
a, b, c, d = A'

@show A
@show oddsratiohat(a, b, c, d)
@show confint_or_pearson_chisq(a, b, c, d)
@show pvalue_or_pearson_chisq(a, b, c, d)
;

# %%
A = [
    11 17969-11
     2 18327-2
]
a, b, c, d = A'

@show A
@show oddsratiohat(a, b, c, d)
@show confint_or_pearson_chisq(a, b, c, d)
@show pvalue_or_pearson_chisq(a, b, c, d)
;

# %%
A = [
    11 17969-11
     3 18327-3
]
a, b, c, d = A'

@show A
@show oddsratiohat(a, b, c, d)
@show confint_or_pearson_chisq(a, b, c, d)#correction=0.5)
@show pvalue_or_pearson_chisq(a, b, c, d)#correction=0.5)
;

# %%
A = [
    11 17969-11
     3 18327-3
]
a, b, c, d = A'

@show A
@show oddsratiohat(a, b, c, d)
@show confint_or_pearson_chisq(a, b, c, d; correction=0.5)
@show pvalue_or_pearson_chisq(a, b, c, d; correction=0.5)
;

# %%
