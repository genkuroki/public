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
default(fmt=:png)

ECDF(A, x) = count(≤(x), A)/length(A)

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
    αs = range(0.002, 1, 1000)
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
safediv(x, y) = x == 0 ? zero(x/y) : x/y

expectval(A) = safediv.(sum(A; dims=2) * sum(A; dims=1), sum(A))

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

# %%
function sim_pearson_chisq_stat(randfunc, param; L=10^4)
    A = randfunc(param)
    r, c = size(A)
    df = (r-1)*(c-1)
    chisq = zeros(L)
    pval = zeros(L)
    for i in 1:L
        A = randfunc(param)
        chisq[i] = pearson_chisq_stat(A)
        pval[i] = ccdf.(Chisq(df), chisq[i])
    end
    (; pval, chisq)
end

# %%
x ⪅ y = x < y || x ≈ y

function pvalue_chisq(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    chi2 = (a+b+c+d)*safediv((a*d - b*c)^2, (a+b)*(c+d)*(a+c)*(b+d))
    ccdf(Chisq(1), chi2)
end

function pvalue_yates(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    N = a+b+c+d
    chi2 = N*safediv(max(0, abs(a*d - b*c) - N/2)^2, (a+b)*(c+d)*(a+c)*(b+d))
    ccdf(Chisq(1), chi2)
end

function pvalue_fisher_central(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    hg = Hypergeometric(a+b, c+d, a+c)
    min(1, 2cdf(hg, a), 2ccdf(hg, a-1))
end

function pvalue_fisher_minlike(A)
    @assert size(A) == (2, 2)
    a, b, c, d = A'
    hg = Hypergeometric(a+b, c+d, a+c)
    pa = pdf(hg, a)
    sum(pdf(hg, x) for x in support(hg) if pdf(hg, x) ⪅ pa)
end

function sim_fisher(randfunc, param; L=10^4)
    pval_central = zeros(L)
    pval_minlike = zeros(L)
    pval_chisq = zeros(L)
    pval_yates = zeros(L)
    Threads.@threads for i in 1:L
        A = randfunc(param)
        pval_central[i] = pvalue_fisher_central(A)
        pval_minlike[i] = pvalue_fisher_minlike(A)
        pval_chisq[i] = pvalue_chisq(A)
        pval_yates[i] = pvalue_yates(A)
    end
    pval_central, pval_minlike, pval_chisq, pval_yates
end

# %%
randpoissons(E) = @. rand(Poisson(E))

@show A = [
    1 2 3
    4 5 6
    7 8 9
]

@show E = expectval(A)

@show [randpoissons(E) for _ in 1:10];

@time (; pval) = sim_pearson_chisq_stat(randpoissons, E; L=10^6)

@show ECDF(pval, 0.05)
@show ECDF(pval, 0.01)
plot_ecdfpval([pval]; linestyles=[:solid])

# %%
function randmultinomial((N, P))
    A = rand(Multinomial(N, vec(P)))
    reshape(A, size(P))
end

@show A = [
    1 2 3
    4 5 6
    7 8 9
]

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N

@show [randmultinomial((N, P)) for _ in 1:10]

@time (; pval) = sim_pearson_chisq_stat(randmultinomial, (N, P); L=10^6)

@show ECDF(pval, 0.05)
@show ECDF(pval, 0.01)
plot_ecdfpval([pval]; linestyles=[:solid])

# %%
function randmultinomials((Ns, Ps))
    r = length(Ns)
    A = similar(Ps, Int)
    for i in 1:r
        rand!(Multinomial(Ns[i], @view(Ps[i,:])), @view(A[i,:]))
    end
    A
end

@show A = [
    1 2 3
    4 5 6
    7 8 9
]

@show E = expectval(A)
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@show [randmultinomials((Ns, Ps)) for _ in 1:10]

@time (; pval) = sim_pearson_chisq_stat(randmultinomials, (Ns, Ps); L=10^6)

@show ECDF(pval, 0.05)
@show ECDF(pval, 0.01)
plot_ecdfpval([pval]; linestyles=[:solid])

# %%
@show A = [
    1 2 3
    4 5 6
    7 8 9
]

@show pearson_chisq_stat(A)
@show pvalue_pearson_chisq(A)
@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time (; pval) = sim_pearson_chisq_stat(randmultinomials, (Ns, Ps); L=10^6)
pval_multinomials = pval
@time (; pval) = sim_pearson_chisq_stat(randmultinomial, (N, P); L=10^6)
pval_multinomial = pval
@time (; pval) = sim_pearson_chisq_stat(randpoissons, E; L=10^6)
pval_poissons = pval

@show ECDF(pval_multinomials, 0.05)
@show ECDF(pval_multinomial, 0.05)
@show ECDF(pval_poissons, 0.05)
@show ECDF(pval_multinomials, 0.01)
@show ECDF(pval_multinomial, 0.01)
@show ECDF(pval_poissons, 0.01)

plot_ecdfpval(
    [pval_multinomials, pval_multinomial, pval_poissons];
    labels=["multinomials", "multinomial", "poissons"])

# %%
@show A = [
     7  4 23
    10  8  8
]

@show pearson_chisq_stat(A)
@show pvalue_pearson_chisq(A)
@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time (; pval) = sim_pearson_chisq_stat(randmultinomials, (Ns, Ps); L=10^6)
pval_multinomials = pval
@time (; pval) = sim_pearson_chisq_stat(randmultinomial, (N, P); L=10^6)
pval_multinomial = pval
@time (; pval) = sim_pearson_chisq_stat(randpoissons, E; L=10^6)
pval_poissons = pval

@show ECDF(pval_multinomials, 0.05)
@show ECDF(pval_multinomial, 0.05)
@show ECDF(pval_poissons, 0.05)
@show ECDF(pval_multinomials, 0.01)
@show ECDF(pval_multinomial, 0.01)
@show ECDF(pval_poissons, 0.01)

plot_ecdfpval(
    [pval_multinomials, pval_multinomial, pval_poissons];
    labels=["multinomials", "multinomial", "poissons"])

# %%
@show A = [
     7  2 23
    10  4  8
]

@show pearson_chisq_stat(A)
@show pvalue_pearson_chisq(A)
@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time (; pval) = sim_pearson_chisq_stat(randmultinomials, (Ns, Ps); L=10^6)
pval_multinomials = pval
@time (; pval) = sim_pearson_chisq_stat(randmultinomial, (N, P); L=10^6)
pval_multinomial = pval
@time (; pval) = sim_pearson_chisq_stat(randpoissons, E; L=10^6)
pval_poissons = pval

@show ECDF(pval_multinomials, 0.05)
@show ECDF(pval_multinomial, 0.05)
@show ECDF(pval_poissons, 0.05)
@show ECDF(pval_multinomials, 0.01)
@show ECDF(pval_multinomial, 0.01)
@show ECDF(pval_poissons, 0.01)

plot_ecdfpval(
    [pval_multinomials, pval_multinomial, pval_poissons];
    labels=["multinomials", "multinomial", "poissons"])

# %%
@show A = [
     7  4 23
     5  4  4
]

@show pearson_chisq_stat(A)
@show pvalue_pearson_chisq(A)
@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time (; pval) = sim_pearson_chisq_stat(randmultinomials, (Ns, Ps); L=10^6)
pval_multinomials = pval
@time (; pval) = sim_pearson_chisq_stat(randmultinomial, (N, P); L=10^6)
pval_multinomial = pval
@time (; pval) = sim_pearson_chisq_stat(randpoissons, E; L=10^6)
pval_poissons = pval

@show ECDF(pval_multinomials, 0.05)
@show ECDF(pval_multinomial, 0.05)
@show ECDF(pval_poissons, 0.05)
@show ECDF(pval_multinomials, 0.01)
@show ECDF(pval_multinomial, 0.01)
@show ECDF(pval_poissons, 0.01)

plot_ecdfpval(
    [pval_multinomials, pval_multinomial, pval_poissons];
    labels=["multinomials", "multinomial", "poissons"])

# %%
@show A = [
    7 2
    2 9
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randpoissons, A; L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 6
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randpoissons, A; L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 6
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randpoissons, E; L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 9
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randmultinomial, (N, A/N); L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 6
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randmultinomial, (N, A/N); L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 6
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randmultinomial, (N, P); L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 9
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randmultinomials, (Ns, A ./ Ns); L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 6
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randmultinomials, (Ns, A ./ Ns); L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
@show A = [
    7 2
    2 6
]

@show pvalue_chisq(A)
@show pvalue_yates(A)
@show pvalue_fisher_central(A)
@show pvalue_fisher_minlike(A)

@show E = expectval(A)
@show N = sum(A)
@show P = E ./ N
@show Ns = sum(A; dims=2)
@show Ps = E ./ Ns

@time pval_fisher_central, pval_fishet_minlike, pval_chisq, pval_yates = sim_fisher(randmultinomials, (Ns, Ps); L=10^6)

@show ECDF(pval_chisq, 0.05)
@show ECDF(pval_yates, 0.05)
@show ECDF(pval_fisher_central, 0.05)
@show ECDF(pval_fishet_minlike, 0.05)
@show ECDF(pval_chisq, 0.01)
@show ECDF(pval_yates, 0.01)
@show ECDF(pval_fisher_central, 0.01)
@show ECDF(pval_fishet_minlike, 0.01)

plot_ecdfpval(
    [pval_chisq, pval_yates, pval_fisher_central, pval_fishet_minlike];
    labels=["chisq", "yates", "fisher central", "fisher minlike"])

# %%
