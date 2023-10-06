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
using HypothesisTests
using StatsPlots
default(fmt=:png, titlefontsize=12)

ECDF(A, x) = count(≤(x), A)/length(A)

# %%
@show distx = Normal(0, 1)
@show disty = Normal(0, 1.5)
@show m = 25
@show n = 15
@show threshold = 0.05

L = 10^6
pval_Ftest = zeros(L)
pval_Student = zeros(L)
pval_Welch = zeros(L)
pval_multiple = zeros(L)

Threads.@threads for i in 1:L
    X = rand(distx, m)
    Y = rand(disty, n)
    pval_Ftest[i] = pvalue(VarianceFTest(X, Y))
    pval_Student[i] = pvalue(EqualVarianceTTest(X, Y))
    pval_Welch[i] = pvalue(UnequalVarianceTTest(X, Y))
    pval_multiple[i] = pval_Ftest[i] < threshold ? pval_Welch[i] : pval_Student[i]
end

@show ECDF(pval_Ftest, threshold)
@show ECDF(pval_multiple, 0.05)
@show ECDF(pval_Student, 0.05)
@show ECDF(pval_Welch, 0.05)

_tick = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
tick = (_tick, string.(_tick))

P = plot(α -> ECDF(pval_Ftest, α), 0.005, 1; label="")
plot!(xscale=:log10, yscale=:log10)
plot!(xtick=tick, ytick=tick)
plot!(xlim=(0.004, 1), ylim=(0.004, 1))
plot!(xguide="α", yguide="probability of P-value ≤ α")
title!("F-test")
plot!(size=(400, 400))

Q = plot()
plot!(α -> ECDF(pval_multiple, α), 0.005, 1; label="multiple")
plot!(α -> ECDF(pval_Student, α), 0.005, 1; label="Student", ls=:dash)
plot!(α -> ECDF(pval_Welch, α), 0.005, 1; label="Welch", ls=:dashdot)
plot!(identity, 0.005, 1; label="", ls=:dot, c=:black, alpha=0.5)
plot!(xscale=:log10, yscale=:log10)
plot!(xtick=tick, ytick=tick)
plot!(xlim=(0.004, 1), ylim=(0.004, 1))
plot!(xguide="α", yguide="probability of P-value ≤ α")
title!("t-tests")
plot!(size=(400, 400))

plot(P, Q; size=(800, 400))
plot!(leftmargin=4Plots.mm)

# %%
@show distx = Normal(0, 1)
@show disty = Normal(0, 1.5)
@show m = 25
@show n = 15
@show threshold = 0.4

L = 10^6
pval_Ftest = zeros(L)
pval_Student = zeros(L)
pval_Welch = zeros(L)
pval_multiple = zeros(L)

Threads.@threads for i in 1:L
    X = rand(distx, m)
    Y = rand(disty, n)
    pval_Ftest[i] = pvalue(VarianceFTest(X, Y))
    pval_Student[i] = pvalue(EqualVarianceTTest(X, Y))
    pval_Welch[i] = pvalue(UnequalVarianceTTest(X, Y))
    pval_multiple[i] = pval_Ftest[i] < threshold ? pval_Welch[i] : pval_Student[i]
end

@show ECDF(pval_Ftest, threshold)
@show ECDF(pval_multiple, 0.05)
@show ECDF(pval_Student, 0.05)
@show ECDF(pval_Welch, 0.05)

_tick = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
tick = (_tick, string.(_tick))

P = plot(α -> ECDF(pval_Ftest, α), 0.005, 1; label="")
plot!(xscale=:log10, yscale=:log10)
plot!(xtick=tick, ytick=tick)
plot!(xlim=(0.004, 1), ylim=(0.004, 1))
plot!(xguide="α", yguide="probability of P-value ≤ α")
title!("F-test")
plot!(size=(400, 400))

Q = plot()
plot!(α -> ECDF(pval_multiple, α), 0.005, 1; label="multiple")
plot!(α -> ECDF(pval_Student, α), 0.005, 1; label="Student", ls=:dash)
plot!(α -> ECDF(pval_Welch, α), 0.005, 1; label="Welch", ls=:dashdot)
plot!(identity, 0.005, 1; label="", ls=:dot, c=:black, alpha=0.5)
plot!(xscale=:log10, yscale=:log10)
plot!(xtick=tick, ytick=tick)
plot!(xlim=(0.004, 1), ylim=(0.004, 1))
plot!(xguide="α", yguide="probability of P-value ≤ α")
title!("t-tests")
plot!(size=(400, 400))

plot(P, Q; size=(800, 400))
plot!(leftmargin=4Plots.mm)

# %%
