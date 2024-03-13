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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using HypothesisTests
using StatsPlots
default(fmt=:png)
ecdf(A, x) = count(≤(x), A)/length(A)

@show distx = Exponential(1)
@show disty = Exponential(2) - 1
m, n = 20, 30
Niters = 10^6
@show m n Niters
pval = [pvalue(UnequalVarianceTTest(rand(distx, m), rand(disty, n))) for _ in 1:Niters]
plot(α -> ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot, c=:gray)
plot!(xguide="α", yguide="probability of P-value ≤ α", size=(400, 400))

# %%
using Distributions
using HypothesisTests
using StatsPlots
default(fmt=:png)
ecdf(A, x) = count(≤(x), A)/length(A)

@show distx = Bernoulli(0.3)
@show disty = Bernoulli(0.3)
m, n = 20, 30
Niters = 10^6
@show m n Niters
pval = [pvalue(UnequalVarianceTTest(rand(distx, m), rand(disty, n))) for _ in 1:Niters]
plot(α -> ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot, c=:gray)
plot!(xguide="α", yguide="probability of P-value ≤ α", size=(400, 400))

# %%
using Distributions
using HypothesisTests
using StatsPlots
default(fmt=:png)
ecdf(A, x) = count(≤(x), A)/length(A)

m, n = 20, 30
@show distx = Binomial(m, 0.3)
@show disty = Binomial(n, 0.3)
Niters = 10^5
@show Niters
pval1 = [(
            a = rand(distx); b = m - a;
            c = rand(disty); d = n - c;
            pvalue(FisherExactTest(a, b, c, d); method=:central)
        ) for _ in 1:Niters]
plot(α -> ecdf(pval1, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot, c=:gray)
plot!(xguide="α", yguide="probability of P-value ≤ α", size=(400, 400))

# %% tags=[]
using Distributions
using HypothesisTests
using StatsPlots
default(fmt=:png)
ecdf(A, x) = count(≤(x), A)/length(A)

m, n = 20, 30
@show distx = Binomial(m, 0.3)
@show disty = Binomial(n, 0.3)
Niters = 10^5
@show Niters
pval2 = [(
            a = rand(distx); b = m - a;
            c = rand(disty); d = n - c;
            pvalue(FisherExactTest(a, b, c, d); method=:minlike)
        ) for _ in 1:Niters]
plot(α -> ecdf(pval2, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot, c=:gray)
plot!(xguide="α", yguide="probability of P-value ≤ α", size=(400, 400))

# %%
using Distributions
using HypothesisTests
using StatsPlots
default(fmt=:png)
ecdf(A, x) = count(≤(x), A)/length(A)

m, n = 20, 30
@show distx = Binomial(m, 0.3)
@show disty = Binomial(n, 0.3)
Niters = 10^5
@show Niters
@time pval3 = [(
            a = rand(distx); b = m - a;
            c = rand(disty); d = n - c;
            pvalue(ChisqTest([a b; c d]))
        ) for _ in 1:Niters]
plot(α -> ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot, c=:gray)
plot!(xguide="α", yguide="probability of P-value ≤ α", size=(400, 400))

# %%
using Distributions
using HypothesisTests
using StaticArrays
using StatsPlots
default(fmt=:png)
ecdf(A, x) = count(≤(x), A)/length(A)

m, n = 20, 30
@show distx = Binomial(m, 0.3)
@show disty = Binomial(n, 0.3)
Niters = 10^5
@show Niters
@time pval3 = [(
            a = rand(distx); b = m - a;
            c = rand(disty); d = n - c;
            pvalue(ChisqTest(@SMatrix [a b; c d]))
        ) for _ in 1:Niters]
plot(α -> ecdf(pval, α), 0, 0.1; label="")
plot!(identity; label="", ls=:dot, c=:gray)
plot!(xguide="α", yguide="probability of P-value ≤ α", size=(400, 400))

# %%
plot(α -> ecdf(pval1, α), 0, 0.1; label="Fisher central")
plot!(α -> ecdf(pval2, α); label="Fisher minlike")
plot!(α -> ecdf(pval3, α); label="Pearson's χ²")
plot!(identity; label="", ls=:dot, c=:gray)
plot!(xguide="α", yguide="probability of P-value ≤ α", size=(400, 400))

# %%
using Distributions
using HypothesisTests
using StaticArrays
using StatsPlots
default(fmt=:png)
ecdf(A, x) = count(≤(x), A)/length(A)

function plot_sim(; m=20, n=30, p=0.3, q=0.3, Niters=10^5, ytick=0:0.01:1)
    @show distx = Binomial(m, 0.3)
    @show disty = Binomial(n, 0.3)
    @show Niters
    pval1 = zeros(Niters)
    pval2 = zeros(Niters)
    pcal3 = zeros(Niters)
    for i in 1:Niters
        a = rand(distx)
        b = m - a
        c = rand(disty)
        d = n - c
        pval1[i] = pvalue(FisherExactTest(a, b, c, d); method=:central)
        pval2[i] = pvalue(FisherExactTest(a, b, c, d); method=:minlike)
        pval3[i] = pvalue(ChisqTest(@SMatrix [a b; c d]))
    end
    plot(α -> ecdf(pval1, α), 0, 0.1; label="Fisher central")
    plot!(α -> ecdf(pval2, α); label="Fisher minlike", ls=:dash)
    plot!(α -> ecdf(pval3, α); label="Pearson's χ²", ls=:dashdot)
    plot!(identity; label="", ls=:dot, c=:gray)
    plot!(; xguide="α", yguide="probability of P-value ≤ α")
    plot!(; xtick=0:0.01:1, ytick)
    plot!(; size=(400, 400))
end

plot_sim(; m=20, n=30, p=0.3, q=0.3)

# %%
