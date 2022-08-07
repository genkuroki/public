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
using StatsPlots
default(fmt=:png, size=(400, 250))

# %%
function logtick(; xlim=(0.03, 500))
    xmin, xmax = xlim
    a = floor(Int, log10(xmin))
    b = ceil(Int, log10(xmax))
    nums =     [1, 2, 3, 4, 5, 6, 7, 8, 9]
    mask = Bool[1, 1, 0, 0, 1, 0, 0, 0, 0]
    
    logtick = foldl(vcat, ([10.0^k*x for x in nums if xmin ≤ 10.0^k*x ≤ xmax] for k in a:b))
    logticklabel_a = foldl(vcat,
        ([mask[i] ? string(round(10.0^k*x; digits=-k)) : ""
                for (i, x) in enumerate(nums) if xmin ≤ 10.0^k*x ≤ xmax]
            for k in a:-1))
    logticklabel_b = foldl(vcat,
        ([mask[i] ? string(10^k*x) : ""
                for (i, x) in enumerate(nums) if xmin ≤ 10.0^k*x ≤ xmax]
            for k in 0:b))
    logticklabel = vcat(logticklabel_a, logticklabel_b)
    (logtick, logticklabel)
end

logtick()

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
mn = [(m, n) for m in 1:20 for n in 1:2m-1 if !is_symmetric_Hypergeometric(m, m, n)]
@show length(mn)
@show mn;

# %%
@show m, n = mn[20]
fnchg = FisherNoncentralHypergeometric(m, m, n, 1.0)
[(k, pdf(fnchg, k)) for k in support(fnchg)]

# %%
function pvalue_less_equal__J(a, b, c, d; OR = 1.0)
    fnchg = FisherNoncentralHypergeometric(a+b, c+d, a+c, OR)
    sum(pdf(fnchg, x) for x in support(fnchg) if pdf(fnchg, x) ≤ pdf(fnchg, a))
end

x ⪅ y = x < y || x ≈ y

function pvalue_less_approx_J(a, b, c, d; OR = 1.0)
    fnchg = FisherNoncentralHypergeometric(a+b, c+d, a+c, OR)
    sum(pdf(fnchg, x) for x in support(fnchg) if pdf(fnchg, x) ⪅ pdf(fnchg, a))
end

function pvalue_fisher_test_R(a, b, c, d)
    rcopy(stats.fisher_test([a b; c d]))[:p_value]
end

# %%
@show m, n = mn[20]
a = 4
b, c, d = m-a, n-a, m+a-n
@show a, b, c, d
@show pvalue_less_equal__J(a, b, c, d)
@show pvalue_less_approx_J(a, b, c, d)
@show pvalue_fisher_test_R(a, b, c, d);

# %%
@show m, n = mn[23]
fnchg = FisherNoncentralHypergeometric(m, m, n, 1.0)
[(k, pdf(fnchg, k)) for k in support(fnchg)]

# %%
@show m, n = mn[23]
a = 7
b, c, d = m-a, n-a, m+a-n
@show a, b, c, d
@show pvalue_less_equal__J(a, b, c, d)
@show pvalue_less_approx_J(a, b, c, d)
@show pvalue_fisher_test_R(a, b, c, d);

# %%
ORhat = (a/b)/(c/d)

# %%
plot(x -> pvalue_less_equal__J(a, b, c, d; OR = x), 0.1, 500; label="")
plot!(x -> pvalue_less_approx_J(a, b, c, d; OR = x), 0.1, 500; label="", ls=:dash)
vline!([ORhat]; label="")
plot!(xscale=:log10, xtick=logtick(xlim=(0.1, 500)), ytick=0:0.1:1)

# %%
xs = 0.99:0.0005:1.01
plot(xs, x -> pvalue_less_equal__J(a, b, c, d; OR = x); label="")
plot!(xs, x -> pvalue_less_approx_J(a, b, c, d; OR = x); label="", ls=:dash)
plot!(ytick=0:0.005:1)

# %%
m, n = 20, 20
@show m, n
fnchg = FisherNoncentralHypergeometric(m, m, n, 1.0)
[(k, pdf(fnchg, k)) for k in support(fnchg)]

# %%
m, n = 20, 20
@show m, n
a = 13
b, c, d = m-a, n-a, m+a-n
@show a, b, c, d
@show pvalue_less_equal__J(a, b, c, d)
@show pvalue_less_approx_J(a, b, c, d)
@show pvalue_fisher_test_R(a, b, c, d);

# %%
ORhat = (a/b)/(c/d)

# %%
plot(x -> pvalue_less_equal__J(a, b, c, d; OR = x), 0.3, 30; label="")
plot!(x -> pvalue_less_approx_J(a, b, c, d; OR = x), 0.3, 30; label="", ls=:dash)
vline!([ORhat]; label="")
plot!(xscale=:log10, xtick=logtick(xlim=(0.1, 30)), ytick=0:0.1:1)

# %%
xs = 0.99:0.0005:1.01
plot(xs, x -> pvalue_less_equal__J(a, b, c, d; OR = x); label="")
plot!(xs, x -> pvalue_less_approx_J(a, b, c, d; OR = x); label="", ls=:dash)
plot!(ytick=0:0.005:1)

# %%
