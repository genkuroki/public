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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using DataFrames
using Distributions
using RCall

# %%
x ⪅ y = x < y || x ≈ y

function pvalue_fisher_yoshida(a, b, c, d)
    hg = Hypergeometric(a+b, c+d, a+c)
    p_lower = cdf(hg, a)
    p_higher = ccdf(hg, a-1)
    p = 0.0
    if p_lower > p_higher
        for j in support(hg)
            p_j = pdf(hg, j)
            p + p_j ⪅ p_higher || break
            p += p_j
        end
        p += p_higher
    else
        for j in reverse(support(hg))
            p_j = pdf(hg, j)
            p + p_j ⪅ p_lower || break
            p += p_j
        end
        p += p_lower
    end
    p
end

function pvalue_fisher(a, b, c, d)
    hg = Hypergeometric(a+b, c+d, a+c)
    p_a = pdf(hg, a)
    sum(pdf(hg, j) for j in support(hg) if pdf(hg, j) ⪅ p_a)
end

function makedf(a, b, c, d)
    @show a, b, c, d
    hg = Hypergeometric(a+b, c+d, a+c)
    j = reverse(support(hg))
    df = DataFrame(
    j = j, 
        var"j=a" = @.(Int(j == a)), 
        var"P(j)≤P(a)" = @.(Int(pdf(hg, j) ⪅ pdf(hg, a))),
        var"P(j)" = @.(pdf(hg, j)),
        var"P(≥j)" = @.(ccdf(hg, j-1)),
        var"P(≤j)" = @.(cdf(hg, j))
    )
end

# %%
A = [
    12 3
     6 8
]

# %%
makedf(A...)

# %%
@rput A
R"""fisher.test(A)$p.value"""

# %%
R"""
a <- 15
b <- 14
n <- 18
x <- 12
x1 <- c(4:6)
x2 <- c(12:15)
 
sum(
    dhyper(x = x1, m = a, n = b, k = n),
    dhyper(x = x2, m = a, n = b, k = n)
    )
"""

# %%
pvalue_fisher(A...)

# %%
pvalue_fisher_yoshida(A...)

# %%
for _ in 1:1000
    A = rand(2:10, 2, 2)
    if !(pvalue_fisher(A...) ≈ pvalue_fisher_yoshida(A...)) && pvalue_fisher(A...) < 0.3
        @show A
        break
    end
end

# %%
A = [
    7  8
    2 10
]

# %%
makedf(A...)

# %%
@rput A
R"""fisher.test(A)$p.value"""

# %%
R"""
a <- 15
b <- 12
n <- 9
x <- 7
x1 <- c(7:9)
x2 <- c(0:2)
 
sum(
    dhyper(x = x1, m = a, n = b, k = n),
    dhyper(x = x2, m = a, n = b, k = n)
    )
"""

# %%
pvalue_fisher(A...)

# %%
pvalue_fisher_yoshida(A...)

# %%
