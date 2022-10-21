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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
ENV["COLUMNS"] = 120

using HypothesisTests

# minlik版のFisher検定の実装では ≤ ではなく, 次の ⪅ を使う必要がある.
x ⪅ y = x < y || x ≈ y

# minlik版Fisher検定のP値の素朴だが(高率の悪いが)正しい実装
using Distributions
function fisher_test(a, b, c, d)
    hg = Hypergeometric(a+b, c+d, a+c)
    pa = pdf(hg, a)
    sum(pdf(hg, i) for i in support(hg) if pdf(hg, i) ⪅ pa)
end

# ≤ を使う正しくない実装
function fisher_test_le(a, b, c, d)
    hg = Hypergeometric(a+b, c+d, a+c)
    pa = pdf(hg, a)
    sum(pdf(hg, i) for i in support(hg) if pdf(hg, i) ≤ pa)
end

# https://qiita.com/WolfMoon/items/530413ce7c8439d18c18
using SpecialFunctions
# 以下も不適切な実装
function fisher_test_wm(a, b, c, d)
    lchoose(n, k) = logfactorial(n) - logfactorial(k) - logfactorial(n - k)
    Stats(i, e, f, g, n) = exp(lchoose(e, i) + lchoose(f, g - i) - lchoose(n, g))
    e, f, g, h, n = a + b, c + d, a + c, b + d, a + b + c + d
    mi = max(0, e + g - n)
    length = min(e, g) - mi
    prob = [Stats(mi + i, e, f, g, n) for i in 0:length]
    #println("Fisher's Exact Test for Count Data (two tailed)")
    #println("p value =", sum(prob[prob .<= Stats(a, e, f, g, n)]))
    sum(prob[prob .<= Stats(a, e, f, g, n)])
end

# 訂正版
function fisher_test_rev(a, b, c, d)
    lchoose(n, k) = logfactorial(n) - logfactorial(k) - logfactorial(n - k)
    Stats(i, e, f, g, n) = exp(lchoose(e, i) + lchoose(f, g - i) - lchoose(n, g))
    e, f, g, h, n = a + b, c + d, a + c, b + d, a + b + c + d
    mi = max(0, e + g - n)
    length = min(e, g) - mi
    prob = [Stats(mi + i, e, f, g, n) for i in 0:length]
    #println("Fisher's Exact Test for Count Data (two tailed)")
    #println("p value =", sum(prob[prob .<= Stats(a, e, f, g, n)]))
    sum(prob[prob .⪅ Stats(a, e, f, g, n)])
end

# %%
s, t, n = 17, 7, 12
hg = Hypergeometric(s, t, n)
supp = support(hg)
Any[
    [(a, s-a, n-a, t-n+a) for a in supp];;
    [fisher_test_wm(a, s-a, n-a, t-n+a) for a in supp];;
    [pvalue(FisherExactTest(a, s-a, n-a, t-n+a), method=:minlike) for a in supp];;
    [fisher_test(a, s-a, n-a, t-n+a) for a in supp];;
    [fisher_test_le(a, s-a, n-a, t-n+a) for a in supp];;
    [fisher_test_wm(a, s-a, n-a, t-n+a) ≈ pvalue(FisherExactTest(a, s-a, n-a, t-n+a), method=:minlike) for a in supp];;
    [fisher_test_wm(a, s-a, n-a, t-n+a) ≈ fisher_test(a, s-a, n-a, t-n+a) for a in supp];;
    [fisher_test(a, s-a, n-a, t-n+a) ≈ fisher_test_le(a, s-a, n-a, t-n+a) for a in supp];;
    [fisher_test_rev(a, s-a, n-a, t-n+a) for a in supp];;
    [fisher_test_rev(a, s-a, n-a, t-n+a) ≈ fisher_test(a, s-a, n-a, t-n+a) for a in supp];;
]

# %%
a, b, c, d = 11, 6, 1, 6
@show fisher_test_wm(a, b, c, d)
@show pvalue(FisherExactTest(a, b, c, d))
@show fisher_test(a, b, c, d)
@show fisher_test_le(a, b, c, d);

# %%
a, b, c, d = 9, 8, 3, 4
@show fisher_test_wm(a, b, c, d)
@show pvalue(FisherExactTest(a, b, c, d))
@show fisher_test(a, b, c, d)
@show fisher_test_le(a, b, c, d);

# %%
s, t, n = 17, 7, 12
hg = Hypergeometric(s, t, n)
supp = support(hg)
[(a, pdf(Hypergeometric(s, t, n), a)) for a in supp]

# %%
using HypothesisTests
using RCall
using SpecialFunctions

# https://qiita.com/WolfMoon/items/530413ce7c8439d18c18
function fisher(a, b, c, d)
    lchoose(n, k) = logfactorial(n) - logfactorial(k) - logfactorial(n - k)
    Stats(i, e, f, g, n) = exp(lchoose(e, i) + lchoose(f, g - i) - lchoose(n, g))
    e, f, g, h, n = a + b, c + d, a + c, b + d, a + b + c + d
    mi = max(0, e + g - n)
    length = min(e, g) - mi
    prob = [Stats(mi + i, e, f, g, n) for i in 0:length]
    println("Fisher's Exact Test for Count Data (two tailed)")
    println("p value =", sum(prob[prob .<= Stats(a, e, f, g, n)]))
end;

# %%
fisher(11, 6, 1, 6)

# %%
pvalue(FisherExactTest(11, 6, 1, 6), method=:minlike)

# %%
R"fisher.test(matrix(c(11, 6, 1, 6), nrow=2))"

# %%
fisher(9, 8, 3, 4)

# %%
pvalue(FisherExactTest(9, 8, 3, 4), method=:minlike)

# %%
R"fisher.test(matrix(c(9, 8, 3, 4), nrow=2))"

# %%
