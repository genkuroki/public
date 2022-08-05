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
using RCall
@rimport stats as stats

# %%
function is_symmetric_Hypergeometric(s, f, n)
    fnchg = FisherNoncentralHypergeometric(s, f, n, 1.0)
    xmin, xmax = extrema(fnchg)
    for i in 0:(xmax-xmin)÷2
        pdf(fnchg, xmin+i) != pdf(fnchg, xmax-i) && return false
    end
    true
end

# %%
mn = [(m, n) for m in 1:10 for n in 1:2m-1 if !is_symmetric_Hypergeometric(m, m, n)]
@show mn;

# %%
@show m, n = mn[20]
fnchg = FisherNoncentralHypergeometric(m, m, n, 1.0)
[(k, pdf(fnchg, k)) for k in support(fnchg)]

# %%
function pvalue_less_equal_(a, b, c, d)
    fnchg = FisherNoncentralHypergeometric(a+b, c+d, a+c, 1.0)
    sum(pdf(fnchg, x) for x in support(fnchg) if pdf(fnchg, x) ≤ pdf(fnchg, a))
end

x ⪅ y = x < y || x ≈ y

function pvalue_less_approx(a, b, c, d)
    fnchg = FisherNoncentralHypergeometric(a+b, c+d, a+c, 1.0)
    sum(pdf(fnchg, x) for x in support(fnchg) if pdf(fnchg, x) ⪅ pdf(fnchg, a))
end

function pvalue_fisher_test(a, b, c, d)
    rcopy(stats.fisher_test([a b; c d]))[:p_value]
end

# %%
@show m, n = mn[20]
a = 4
b, c, d = m-a, n-a, m+a-n
@show a, b, c, d
@show pvalue_less_equal_(a, b, c, d)
@show pvalue_less_approx(a, b, c, d)
@show pvalue_fisher_test(a, b, c, d);

# %%
