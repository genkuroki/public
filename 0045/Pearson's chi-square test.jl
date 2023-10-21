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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using Random
using StatsBase: ecdf
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6)

ECDF(A, x) = count(≤(x), A)/length(A)
safediv(x, y) = x == 0 ? zero(x/y) : x/y
x ⪅ y = x < y || x ≈ y

# %%
function plot_ecdfpval(pvals;
        n = length(pvals), 
        labels = fill("", n),
        linestyles = fill(:auto, n),
        size = (400, 400),
        legend = :bottomright,
        kwargs...
    )
    _tick = Any[0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    xtick = ytick = (float.(_tick), string.(_tick))
    xlim = ylim = (0.0015, 1.1)
    αs = exp.(range(log(0.002), log(1), 1000))
    P = plot()
    for (pval, label, ls) in zip(pvals, labels, linestyles)
        _ecdf_pval = ecdf(pval)
        plot!(αs, α -> _ecdf_pval(α); label, ls, kwargs...)
    end
    plot!(αs, x->x; label="", ls=:dot, alpha=0.5, c=:black)
    plot!(αs, x->0.8x; label="", ls=:dot, alpha=0.3, c=:black)
    plot!(αs, x->1.2x; label="", ls=:dot, alpha=0.3, c=:black)
    plot!(; xscale=:log10, yscale=:log10, xtick, ytick, xlim, ylim)
    plot!(; xguide="α", yguide="probability of P-value ≤ α")
    plot!(; size, legend)
end

# %%
pval = rand(10^6)
plot_ecdfpval([pval]; labels=["sample of uniform dist."])

# %%
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

function pvalue_minlike(dist::DiscreteUnivariateDistribution, x)
    Px = pdf(dist, x)
    Px == 0 && return Px
    Px == 1 && return Px
    m = mode(dist)
    Px ≈ pdf(dist, m) && return one(Px)
    if x < m
        y = _search_boundary(_pdf_le, 2m - x, 1, (dist, Px))
        cdf(dist, x) + ccdf(dist, y-1)
    else # x > m
        y = _search_boundary(_pdf_le, 2m - x, -1, (dist, Px))
        cdf(dist, y) + ccdf(dist, x-1)
    end
end

function pvalue_score_bin(k, n, p)
    p̂ = k/n
    sehat = √(p*(1-p)/n)
    z = safediv(p̂ - p, sehat)
    2ccdf(Normal(), abs(z))
end

function pvalue_wald_bin(k, n, p)
    p̂ = k/n
    sehat = √(p̂*(1-p̂)/n)
    z = safediv(p̂ - p, sehat)
    2ccdf(Normal(), abs(z))
end

pvalue_minlike_bin(k, n, p) = pvalue_minlike(Binomial(n, p), k)

function pvalue_central_bin(k, n, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

# %%
function plot_bin(; n = 50, p = 0.2, L = 10^6)
    @show bin = Binomial(n, p)
    K = rand(bin, L)
    pval_score = pvalue_score_bin.(K, n, p)
    pval_wald = pvalue_wald_bin.(K, n, p)
    pval_minlike = pvalue_minlike_bin.(K, n, p)
    pval_central = pvalue_central_bin.(K, n, p)

    X = [
        (pval_score, "score test", 1),
        (pval_wald, "Wald test", 2),
        (pval_minlike, "exact test (minlike)", 3),
        (pval_central, "exact test (central)", 4),
    ]
    PP = []
    for (pval, name, c) in X
        P = plot_ecdfpval([pval]; linestyles=[:solid], c)
        title!("binomial $name")
        push!(PP, P)
    end

    plot(PP...; size=(800, 800), layout=(2, 2))
end

# %%
plot_bin(n = 20, p = 0.2)

# %%
plot_bin(n = 80, p = 0.05)

# %%
randpoissons(E) = @. rand(Poisson(E))

function randmultinomial((N, P))
    A = rand(Multinomial(N, vec(P)))
    reshape(A, size(P))
end

function randmultinomials((Ns, Ps))
    r = length(Ns)
    A = similar(Ps, Int)
    for i in 1:r
        rand!(Multinomial(Ns[i], @view(Ps[i,:])), @view(A[i,:]))
    end
    A
end

# %%
expectval(A) = sum(A; dims=2) * sum(A; dims=1) / sum(A)

function pvalue_pearson_chisq_2x2(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    chi2 = (a+b+c+d)*safediv((a*d - b*c)^2, (a+b)*(c+d)*(a+c)*(b+d))
    ccdf(Chisq(1), chi2)
end

function pvalue_yates_chisq_2x2(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    N = a+b+c+d
    chi2 = N*safediv(max(0, abs(a*d - b*c) - N/2)^2, (a+b)*(c+d)*(a+c)*(b+d))
    ccdf(Chisq(1), chi2)
end

function pvalue_fisher_central_2x2(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    hg = Hypergeometric(a+b, c+d, a+c)
    min(1, 2cdf(hg, a), 2ccdf(hg, a-1))
end

function pvalue_fisher_minlike_2x2_slow(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    hg = Hypergeometric(a+b, c+d, a+c)
    pa = pdf(hg, a)
    sum(pdf(hg, x) for x in support(hg) if pdf(hg, x) ⪅ pa)
end

function pvalue_fisher_minlike_2x2(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    hg = Hypergeometric(a+b, c+d, a+c)
    pvalue_minlike(hg, a)
end

function sim_2x2(randfunc, param; L=10^6)
    pval_fisher_central = zeros(L)
    pval_fisher_minlike = zeros(L)
    pval_pearson_chisq = zeros(L)
    pval_yates_chisq = zeros(L)
    Threads.@threads for i in 1:L
        A = randfunc(param)
        pval_fisher_central[i] = pvalue_fisher_central_2x2(A)
        pval_fisher_minlike[i] = pvalue_fisher_minlike_2x2(A)
        pval_pearson_chisq[i] = pvalue_pearson_chisq_2x2(A)
        pval_yates_chisq[i] = pvalue_yates_chisq_2x2(A)
    end
    pval_fisher_central, pval_fisher_minlike, pval_pearson_chisq, pval_yates_chisq
end

function plot_2x2(A; L=10^6)
    println("A =")
    show(stdout, MIME("text/plain"), A)
    println()
    println()
    
    @show pvalue_pearson_chisq_2x2(A)
    @show pvalue_yates_chisq_2x2(A)
    @show pvalue_fisher_central_2x2(A)
    @show pvalue_fisher_minlike_2x2(A)
    println()

    E = expectval(A)
    println("E = expectval(A) =")
    show(stdout, MIME("text/plain"), E)
    println()
    println()
    @show N = sum(A)
    @show Ns = sum(A; dims=2)
    println()
    
    randfuncs = (
        randpoissons, randpoissons,
        randmultinomial, randmultinomial,
        randmultinomials, randmultinomials,
    )
    params = (
        E, A,
        (N, E / N), (N, A / N),
        (Ns, E ./ Ns), (Ns, A ./ Ns),
    )
    names = (
        "4 Poissons under the null", "4 Poissons with expectation A", 
        "Multinomial under the null", "Multinomial with expectation A", 
        "2 binomials under the null", "2 binomials with expectation A", 
    )

    PP = []
    for (randfunc, param, name) in zip(randfuncs, params, names)
        (
            pval_fisher_central,
            pval_fisher_minlike,
            pval_pearson_chisq,
            pval_yates_chisq
        ) = sim_2x2(randfunc, param; L)
        println("-"^20, " ", "Model: $name")
        @show ECDF(pval_pearson_chisq, 0.05)
        @show ECDF(pval_yates_chisq, 0.05)
        @show ECDF(pval_fisher_central, 0.05)
        @show ECDF(pval_fisher_minlike, 0.05)
        println()
        @show ECDF(pval_pearson_chisq, 0.01)
        @show ECDF(pval_yates_chisq, 0.01)
        @show ECDF(pval_fisher_central, 0.01)
        @show ECDF(pval_fisher_minlike, 0.01)
        println()
        P = plot_ecdfpval(
            [pval_pearson_chisq, pval_yates_chisq, pval_fisher_central, pval_fisher_minlike];
            labels=["Pearson χ²", "Yates χ²", "Fisher central", "Fisher minlike"])
        title!("Model: $name")
        push!(PP, P)
    end
    
    plot(PP...; size=(820, 1200), layout=(3, 2))
    plot!(titlefontsize=10, tickfontsize=6, leftmargin=6Plots.mm)
end

# %%
A = [
    10000 10000
    10000 10400
]

p1 = pvalue_fisher_minlike_2x2_slow(A)
@time p1 = pvalue_fisher_minlike_2x2_slow(A)
@time p1 = pvalue_fisher_minlike_2x2_slow(A)
@time p1 = pvalue_fisher_minlike_2x2_slow(A)
p2 = pvalue_fisher_minlike_2x2(A)
@time p2 = pvalue_fisher_minlike_2x2(A)
@time p2 = pvalue_fisher_minlike_2x2(A)
@time p2 = pvalue_fisher_minlike_2x2(A)
@show p1 p2
p1 ≈ p2

# %%
plot_2x2([3 1; 1 3])

# %%
plot_2x2([6 2; 2 7])

# %%
plot_2x2([9 2; 2 7])

# %%
plot_2x2([19 9; 4 9])

# %%
function pearson_chisq_stat(A)
    E = expectval(A)
    sum(safediv((a - e)^2, e) for (a, e) in zip(A, E))
end

function pvalue_pearson_chisq(A)
    r, c = size(A)
    df = (r-1)*(c-1)
    chi2 = pearson_chisq_stat(A)
    ccdf(Chisq(df), chi2)
end

function sim_pearson_chisq_test(randfunc, param; L=10^6)
    pval = zeros(L)
    for i in 1:L
        A = randfunc(param)
        pval[i] = pvalue_pearson_chisq(A)
    end
    pval
end

function plot_pearson_chisq(A; L=10^6)
    println("A =")
    show(stdout, MIME("text/plain"), A)
    println()
    println()
    
    @show pvalue_pearson_chisq(A)
    println()

    E = expectval(A)
    println("E = expectval(A) =")
    show(stdout, MIME("text/plain"), E)
    println()
    println()
    @show N = sum(A)
    @show Ns = sum(A; dims=2)
    println()
    
    randfuncs = (
        randpoissons, randpoissons,
        randmultinomial, randmultinomial,
        randmultinomials, randmultinomials,
    )
    params = (
        E, A,
        (N, E / N), (N, A / N),
        (Ns, E ./ Ns), (Ns, A ./ Ns),
    )
    names = (
        "Poissons under the null", "Poissons with expectation A", 
        "Multinomial under the null", "Multinomial with expectation A", 
        "Multinomials under the null", "Multinomials with expectation A", 
    )

    PP = []
    for (randfunc, param, name) in zip(randfuncs, params, names)
        pval = sim_pearson_chisq_test(randfunc, param; L)
        println("-"^20, " ", "Model: $name")
        @show ECDF(pval, 0.05)
        @show ECDF(pval, 0.01)
        P = plot_ecdfpval([pval]; linestyles=[:solid])
        title!("Model: $name")
        push!(PP, P)
    end
    println()
    
    plot(PP...; size=(820, 1200), layout=(3, 2))
    plot!(titlefontsize=10, leftmargin=6Plots.mm)
end

# %%
A = [
     7  4 23
    10  8  8
]
plot_pearson_chisq(A)

# %%
A = [
     7  2 23
    10  4  8
]
plot_pearson_chisq(A)

# %%
A = [
     7  4 23
     5  4  4
]
plot_pearson_chisq(A)

# %%
