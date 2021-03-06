# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
#  http://www.snap-tck.com/room04/c01/stat/stat01/stat0106.html#note03

# %%
using Logging
disable_logging(Logging.Warn)

using Distributions
using Roots
using Memoization
using Plots
using RCall
@rlibrary stats

x ⪅ y = x < y || x ≈ y

function pvalue(d, k, ::Val{:ts}) # :ts stands for "two-sided"
    sum(pdf(d, j) for j in support(d) if pdf(d, j) ⪅ pdf(d, k); init = 0.0)
end

function pvalue(d, k, ::Val{:dos}) # :dos stands for "doubled one-sided"
    min(1, 2cdf(d, k), 2ccdf(d, k-1))
end

@memoize function pvalue(a, b, c, d, ::Val{:fisher}; ω = 1.0)
    (iszero(a+c) || iszero(b+d) || iszero(a+b) || iszero(c+d)) && return 1.0
    pvalue(FisherNoncentralHypergeometric(a+b, c+d, a+c, ω), a, Val(:ts))
end

@memoize function pvalue(a, b, c, d, ::Val{:fisher_dos}; ω = 1.0)
    (iszero(a+c) || iszero(b+d) || iszero(a+b) || iszero(c+d)) && return 1.0
    pvalue(FisherNoncentralHypergeometric(a+b, c+d, a+c, ω), a, Val(:dos))
end

function delta(a, b, c, d, ω=1.0)
    A = 1 - ω
    B = a + d + ω*(b + c)
    C = a*d - ω*b*c
    2*C/(-B - √(B^2 - 4*A*C))
end

_odds_ratio(a, b, c, d) = (ad = a*d; iszero(ad) ? ad : ad/(b*c))
function odds_ratio(a, b, c, d, ω=1.0)
    δ = delta(a, b, c, d, ω)
    _odds_ratio(a + δ, b - δ, c - δ, d + δ)
end

_chisq(a, b, c, d) = (a+b+c+d)*(a*d - b*c)^2/((a+b)*(c+d)*(a+c)*(b+d))
function chisq(a, b, c, d, ω=1.0)
    δ = delta(a, b, c, d, ω)
    iszero(δ) ? δ : δ^2 * (1/(a + δ) + 1/(b - δ) + 1/(c - δ) + 1/(d + δ))
end

function chisq_yates(a, b, c, d, ω=1.0)
    δ = delta(a, b, c, d, ω)
    m = max(0, abs(δ)-1/2)^2
    iszero(m) ? m : m * (1/(a + δ) + 1/(b - δ) + 1/(c - δ) + 1/(d + δ))
end

@memoize function pvalue(a, b, c, d, ::Val{:chisq}; ω = 1.0)
    ccdf(Chisq(1), chisq(a, b, c, d, ω))
end

@memoize function pvalue(a, b, c, d, ::Val{:chisq_yates}; ω = 1.0)
    ccdf(Chisq(1), chisq_yates(a, b, c, d, ω))
end

@memoize function confidence_interval(a, b, c, d, alg; α = 0.05)
    CI = exp.(find_zeros(t -> pvalue(a, b, c, d, alg; ω = exp(t)) - α, -20, 20))
    or = _odds_ratio(a, b, c, d)
    iszero(or) && return (0.0, CI[1])
    isinf(or)  && return (CI[1], Inf)
    isone(length(CI)) ? CI[1] < 1 ? (0.0, CI[1]) : (CI[1], Inf) : (CI[1], CI[end])
end

# %%
println("="^30 * " Test data:\n")

A = [16 4; 4 6]
@show A
println()

println("="^30 * " Calculated by R:\n")
@show rcopy(fisher_test(A))[:p_value]
@show rcopy(fisher_test(A))[:conf_int]
println()
@show rcopy(chisq_test(A, correct=false))[:p_value]
@show rcopy(chisq_test(A))[:p_value]

flush(stdout)
println()

println("="^30 * " Calculated by Julia:\n")
@show pvalue(A..., Val(:fisher))
@show confidence_interval(A..., Val(:fisher))
@show pvalue(A..., Val(:fisher_dos))
@show confidence_interval(A..., Val(:fisher_dos))
println()
@show pvalue(A..., Val(:chisq))
@show confidence_interval(A..., Val(:chisq))
@show pvalue(A..., Val(:chisq_yates))
@show confidence_interval(A..., Val(:chisq_yates))
;

# %%
fisher_test(A)

# %%
chisq_test(A, correct = false)

# %%
chisq_test(A)

# %%
function dnhyper(a, b, c, d, ncp)
    nchg = FisherNoncentralHypergeometric(a+b, c+d, a+c, ncp)
    supp = support(nchg)
    pdf.(nchg, supp)
end

# https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/fisher.test.R#L132
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html
function dnhyper_R(a, b, c, d, ncp)
    hg = Hypergeometric(a+b, c+d, a+c)
    supp = support(hg)
    logdc = logpdf.(hg, supp)
    d = @. logdc + log(ncp) * supp
    m = maximum(d)
    d = @. exp(d - m)
    d ./ sum(d)
end

f(ncp) = dnhyper(A..., ncp)
g(ncp) = dnhyper_R(A..., ncp)
@show round.(f(0.5) ./ g(0.5); digits=5)
@show round.(f(1.0) ./ g(1.0); digits=5)
@show round.(f(2.0) ./ g(2.0); digits=5);

# %%
# https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/fisher.test.R#L208
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/uniroot
# https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/nlm.R#L60
# `tol = .Machine$double.eps^0.25`

@show rcopy(R".Machine$double.eps^0.25")
@show d = eps()^0.25
@show CI_R = rcopy(fisher_test(A))[:conf_int]
@show CI = confidence_interval(A..., Val(:fisher_dos))
@show log.(CI_R ./ CI);

# %%
ecdf(A, x) = count(≤(x), A)/length(A)

function multinomial_null(a, b, c, d)
    n = a + b + c + d
    pa = (a+b)*(a+c)/n^2
    pb = (a+b)*(b+d)/n^2
    pc = (c+d)*(a+c)/n^2
    pd = (c+d)*(b+d)/n^2
    Multinomial(n, [pa, pb, pc, pd])
end

null = multinomial_null(A...)
L = 10^5
X = rand(null, L)
@time pvalue_fisher = vec(mapslices(A -> pvalue(A..., Val(:fisher)), X; dims=1))
@time pvalue_fisher_dos = vec(mapslices(A -> pvalue(A..., Val(:fisher_dos)), X; dims=1))
@time pvalue_chisq = vec(mapslices(A -> pvalue(A..., Val(:chisq)), X; dims=1))
@time pvalue_chisq_yates = vec(mapslices(A -> pvalue(A..., Val(:chisq_yates)), X; dims=1));

# %%
p = range(0, 1; length=1001)
P = plot()
plot!(p, p -> ecdf(pvalue_fisher, p); label="fisher")
plot!(p, p -> ecdf(pvalue_fisher_dos, p); label="fisher_dos", ls=:dash)
plot!(p, p -> ecdf(pvalue_chisq, p); label="chisq")
plot!(p, p -> ecdf(pvalue_chisq_yates, p); label="chisq_yates", ls=:dashdot)
plot!([0, maximum(p)], [0, maximum(p)]; label="", color=:black, ls=:dot)
plot!(; xlim=extrema(p), ylim=extrema(p))
plot!(; size=(400, 400), legend=:topleft)
plot!(; xlabel="significance level α", ylabel="probability of the first class errors")
plot!(; xtick=range(extrema(p)...; length=11), ytick=range(extrema(p)...; length=11), tickfontsize=7, guidefontsize=9)
plot!(; title="null = Multinomial(n=$(ntrials(null)), p=$(round.(probs(null); digits=2)))", titlefontsize=9)

p = range(0, 0.1; length=101)
Q = plot()
plot!(p, p -> ecdf(pvalue_fisher, p); label="fisher")
plot!(p, p -> ecdf(pvalue_fisher_dos, p); label="fisher_dos", ls=:dash)
plot!(p, p -> ecdf(pvalue_chisq, p); label="chisq")
plot!(p, p -> ecdf(pvalue_chisq_yates, p); label="chisq_yates", ls=:dashdot)
plot!([0, maximum(p)], [0, maximum(p)]; label="", color=:black, ls=:dot)
plot!(; xlim=extrema(p), ylim=extrema(p))
plot!(; size=(400, 400), legend=:topleft)
plot!(; xlabel="significance level α", ylabel="probability of the first class errors")
plot!(; xtick=range(extrema(p)...; length=11), ytick=range(extrema(p)...; length=11), tickfontsize=7, guidefontsize=9)
plot!(; title="null = Multinomial(n=$(ntrials(null)), p=$(round.(probs(null); digits=2)))", titlefontsize=9)

plot(P, Q; size=(800, 400), left_margin=5Plots.mm, bottom_margin=5Plots.mm)

# %%
logω = range(-0.5, 2.0; length=500)
P = plot()
plot!(logω, t -> pvalue(A..., Val(:fisher); ω=10^t); label="fisher")
plot!(logω, t -> pvalue(A..., Val(:fisher_dos); ω=10^t); label="fisher_dos", ls=:dash)
plot!(logω, t -> pvalue(A..., Val(:chisq); ω=10^t); label="chisq")
plot!(logω, t -> pvalue(A..., Val(:chisq_yates); ω=10^t); label="chisq_yates", ls=:dashdot)
plot!(logω, 0.05ones(length(logω)); label="α=0.05", color=:red)
#plot!(; legend=:topleft)
x = [0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
plot!(; xtick = (log10.(x), string.(x)), ytick=0:0.1:1)
plot!(; title="p-value functions for A = $A", titlefontsize=11)

Q = deepcopy(P)
plot!(Q; ylim=(0, 0.15), ytick=0:0.05:1)

plot(P, Q; size=(800, 300))

# %% [markdown]
# https://twitter.com/myuuuuun/status/1428579692219490304?s=21

# %%
println("="^30 * " Test data:\n")

a, b, c, d = 83207-1557, 1557, 317-1, 1
A = [a b; c d]
@show A
println()

println("="^30 * " Calculated by R:\n")
@show rcopy(fisher_test(A))[:p_value]
@show rcopy(fisher_test(A))[:conf_int]
println()
@show rcopy(chisq_test(A, correct=false))[:p_value]
@show rcopy(chisq_test(A))[:p_value]

sleep(0.1)
println()

println("="^30 * " Calculated by Julia:\n")
@show pvalue(A..., Val(:fisher))
@show confidence_interval(A..., Val(:fisher))
@show pvalue(A..., Val(:fisher_dos))
@show confidence_interval(A..., Val(:fisher_dos))
println()
@show pvalue(A..., Val(:chisq))
@show confidence_interval(A..., Val(:chisq))
@show pvalue(A..., Val(:chisq_yates))
@show confidence_interval(A..., Val(:chisq_yates))
;

# %%
a, b, c, d = 83207-1557, 1557, 317-1, 1
A = [a b; c d]

null = multinomial_null(A...)
L = 10^5
X = rand(null, L)
@time pvalue_fisher = vec(mapslices(A -> pvalue(A..., Val(:fisher)), X; dims=1))
@time pvalue_fisher_dos = vec(mapslices(A -> pvalue(A..., Val(:fisher_dos)), X; dims=1))
@time pvalue_chisq = vec(mapslices(A -> pvalue(A..., Val(:chisq)), X; dims=1))
@time pvalue_chisq_yates = vec(mapslices(A -> pvalue(A..., Val(:chisq_yates)), X; dims=1));


# %%
p = range(0, 1; length=1001)
P = plot()
plot!(p, p -> ecdf(pvalue_fisher, p); label="fisher")
plot!(p, p -> ecdf(pvalue_fisher_dos, p); label="fisher_dos", ls=:dash)
plot!(p, p -> ecdf(pvalue_chisq, p); label="chisq")
plot!(p, p -> ecdf(pvalue_chisq_yates, p); label="chisq_yates", ls=:dashdot)
plot!([0, maximum(p)], [0, maximum(p)]; label="", color=:black, ls=:dot)
plot!(; xlim=extrema(p), ylim=extrema(p))
plot!(; size=(400, 400), legend=:topleft)
plot!(; xlabel="significance level α", ylabel="probability of the first class errors")
plot!(; xtick=range(extrema(p)...; length=11), ytick=range(extrema(p)...; length=11), tickfontsize=7, guidefontsize=9)
plot!(; title="null = Multinomial(n=$(ntrials(null)), p=$(round.(probs(null); digits=4)))", titlefontsize=7)

p = range(0, 0.1; length=101)
Q = plot()
plot!(p, p -> ecdf(pvalue_fisher, p); label="fisher")
plot!(p, p -> ecdf(pvalue_fisher_dos, p); label="fisher_dos", ls=:dash)
plot!(p, p -> ecdf(pvalue_chisq, p); label="chisq")
plot!(p, p -> ecdf(pvalue_chisq_yates, p); label="chisq_yates", ls=:dashdot)
plot!([0, maximum(p)], [0, maximum(p)]; label="", color=:black, ls=:dot)
plot!(; xlim=extrema(p), ylim=extrema(p))
plot!(; size=(400, 400), legend=:topleft)
plot!(; xlabel="significance level α", ylabel="probability of the first class errors")
plot!(; xtick=range(extrema(p)...; length=11), ytick=range(extrema(p)...; length=11), tickfontsize=7, guidefontsize=9)
plot!(; title="null = Multinomial(n=$(ntrials(null)), p=$(round.(probs(null); digits=4)))", titlefontsize=7)

plot(P, Q; size=(800, 400), left_margin=5Plots.mm, bottom_margin=5Plots.mm)

# %%
logω = range(-2.0, 0.5; length=500)
P = plot()
plot!(logω, t -> pvalue(A..., Val(:fisher); ω=10^t); label="fisher")
plot!(logω, t -> pvalue(A..., Val(:fisher_dos); ω=10^t); label="fisher_dos", ls=:dash)
plot!(logω, t -> pvalue(A..., Val(:chisq); ω=10^t); label="chisq")
plot!(logω, t -> pvalue(A..., Val(:chisq_yates); ω=10^t); label="chisq_yates", ls=:dashdot)
plot!(logω, 0.05ones(length(logω)); label="α=0.05", color=:red)
#plot!(; legend=:topleft)
x = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
plot!(; xtick = (log10.(x), string.(x)), ytick=0:0.1:1)
plot!(; title="p-value functions for A = $A", titlefontsize=9)

logω = range(-0.5, 0.5; length=500)
Q = plot()
plot!(logω, t -> pvalue(A..., Val(:fisher); ω=10^t); label="fisher")
plot!(logω, t -> pvalue(A..., Val(:fisher_dos); ω=10^t); label="fisher_dos", ls=:dash)
plot!(logω, t -> pvalue(A..., Val(:chisq); ω=10^t); label="chisq")
plot!(logω, t -> pvalue(A..., Val(:chisq_yates); ω=10^t); label="chisq_yates", ls=:dashdot)
plot!(logω, 0.05ones(length(logω)); label="α=0.05", color=:red)
#plot!(; legend=:topleft)
x = [0.3, 1.0, 3.0]
plot!(; xtick = (log10.(x), string.(x)), ylim=(0, 0.1), ytick=0:0.01:1)
plot!(; title="p-value functions for A = $A", titlefontsize=9)

plot(P, Q; size=(800, 300))

# %% [markdown]
# * https://www3.nhk.or.jp/kansai-news/20210921/2000051589.html
# * https://twitter.com/holyeightturtle/status/1440708315030306824

# %%
println("="^30 * " Test data:\n")

A = [
    14-0   0
    833-78 78
]
@show A
println()

println("="^30 * " Calculated by R:\n")
@show rcopy(fisher_test(A))[:p_value]
@show rcopy(fisher_test(A))[:conf_int]
println()
@show rcopy(chisq_test(A, correct=false))[:p_value]
@show rcopy(chisq_test(A))[:p_value]

flush(stdout)
println()

println("="^30 * " Calculated by Julia:\n")
@show pvalue(A..., Val(:fisher))
@show confidence_interval(A..., Val(:fisher))
@show pvalue(A..., Val(:fisher_dos))
@show confidence_interval(A..., Val(:fisher_dos))
println()
@show pvalue(A..., Val(:chisq))
@show confidence_interval(A..., Val(:chisq))
@show pvalue(A..., Val(:chisq_yates))
@show confidence_interval(A..., Val(:chisq_yates))
;

# %%
A = [
    14-0   0
    833-78 78
]

null = multinomial_null(A...)
L = 10^5
X = rand(null, L)
@time pvalue_fisher = vec(mapslices(A -> pvalue(A..., Val(:fisher)), X; dims=1))
@time pvalue_fisher_dos = vec(mapslices(A -> pvalue(A..., Val(:fisher_dos)), X; dims=1))
@time pvalue_chisq = vec(mapslices(A -> pvalue(A..., Val(:chisq)), X; dims=1))
@time pvalue_chisq_yates = vec(mapslices(A -> pvalue(A..., Val(:chisq_yates)), X; dims=1));

# %%
p = range(0, 1; length=1001)
P = plot()
plot!(p, p -> ecdf(pvalue_fisher, p); label="fisher")
plot!(p, p -> ecdf(pvalue_fisher_dos, p); label="fisher_dos", ls=:dash)
plot!(p, p -> ecdf(pvalue_chisq, p); label="chisq")
plot!(p, p -> ecdf(pvalue_chisq_yates, p); label="chisq_yates", ls=:dashdot)
plot!([0, maximum(p)], [0, maximum(p)]; label="", color=:black, ls=:dot)
plot!(; xlim=extrema(p), ylim=extrema(p))
plot!(; size=(400, 400), legend=:topleft)
plot!(; xlabel="significance level α", ylabel="probability of the first class errors")
plot!(; xtick=range(extrema(p)...; length=11), ytick=range(extrema(p)...; length=11), tickfontsize=7, guidefontsize=9)
plot!(; title="null = Multinomial(n=$(ntrials(null)), p=$(round.(probs(null); digits=4)))", titlefontsize=7)

p = range(0, 0.1; length=101)
Q = plot()
plot!(p, p -> ecdf(pvalue_fisher, p); label="fisher")
plot!(p, p -> ecdf(pvalue_fisher_dos, p); label="fisher_dos", ls=:dash)
plot!(p, p -> ecdf(pvalue_chisq, p); label="chisq")
plot!(p, p -> ecdf(pvalue_chisq_yates, p); label="chisq_yates", ls=:dashdot)
plot!([0, maximum(p)], [0, maximum(p)]; label="", color=:black, ls=:dot)
plot!(; xlim=extrema(p), ylim=extrema(p))
plot!(; size=(400, 400), legend=:topleft)
plot!(; xlabel="significance level α", ylabel="probability of the first class errors")
plot!(; xtick=range(extrema(p)...; length=11), ytick=range(extrema(p)...; length=11), tickfontsize=7, guidefontsize=9)
plot!(; title="null = Multinomial(n=$(ntrials(null)), p=$(round.(probs(null); digits=4)))", titlefontsize=7)

plot(P, Q; size=(800, 400), left_margin=5Plots.mm, bottom_margin=5Plots.mm)

# %%
logω = range(-2.0, 0.5; length=500)
P = plot()
plot!(logω, t -> pvalue(A..., Val(:fisher); ω=10^t); label="fisher")
plot!(logω, t -> pvalue(A..., Val(:fisher_dos); ω=10^t); label="fisher_dos", ls=:dash)
plot!(logω, t -> pvalue(A..., Val(:chisq); ω=10^t); label="chisq")
plot!(logω, t -> pvalue(A..., Val(:chisq_yates); ω=10^t); label="chisq_yates", ls=:dashdot)
plot!(logω, 0.05ones(length(logω)); label="α=0.05", color=:red)
plot!(; legend=:topleft)
x = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]
plot!(; xtick = (log10.(x), string.(x)), ytick=0:0.1:1)
plot!(; title="p-value functions for A = $A", titlefontsize=9)

logω = range(-1.0, 0.1; length=500)
Q = plot()
plot!(logω, t -> pvalue(A..., Val(:fisher); ω=10^t); label="fisher")
plot!(logω, t -> pvalue(A..., Val(:fisher_dos); ω=10^t); label="fisher_dos", ls=:dash)
plot!(logω, t -> pvalue(A..., Val(:chisq); ω=10^t); label="chisq")
plot!(logω, t -> pvalue(A..., Val(:chisq_yates); ω=10^t); label="chisq_yates", ls=:dashdot)
plot!(logω, 0.05ones(length(logω)); label="α=0.05", color=:red)
plot!(; legend=:topleft)
x = [0.3, 1.0, 3.0]
plot!(; xtick = (log10.(x), string.(x)), ylim=(0, 0.1), ytick=0:0.01:1)
plot!(; title="p-value functions for A = $A", titlefontsize=9)

plot(P, Q; size=(800, 300))

# %%
