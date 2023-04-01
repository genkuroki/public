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
using Distributions
E(f, dist) = sum(x -> f(x)*pdf(dist, x), support(dist))
mu(k, dist; m=mean(dist)) = E(x ->(x-m)^k, dist)
using HypothesisTests
using Random
using StatsBase
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6,
    size=(400, 250))

# %%
function sim(distX, distY, m, n; L=10^5)
    pval = Vector{Float64}(undef, L)
    tmpX = [Vector{Float64}(undef, m) for _ in 1:Threads.nthreads()]
    tmpY = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(distX, tmpX[tid])
        Y = rand!(distY, tmpY[tid])
        pval[i] = pvalue(MannWhitneyUTest(X, Y))
    end
    _ecdf_pval = ecdf(pval)
    ecdf_pval(x) = _ecdf_pval(x)
    ecdf_pval
end

# %%
dist7(a, b, c) = Categorical(a, b, c, 1-2(a+b+c), c, b ,a)
@show distX = dist7(1/16, 0, 7/16)
s2, m4, m6 = @show var(distX), mu(4, distX), mu(6, distX)
println()

P = bar(distX; label="", ylim=(0, 0.62), c=1)
title!("distribution of X")

# %%
t = 1/30
a = t
b = (m4 - s2 - 2(3^4-3^2)*a)/(2(2^4 - 2^2))
c = s2/2 - (3^2*a + 2^2*b)
@show distY = dist7(a, b, c)
@show var(distY), mu(4, distY), mu(6, distY)
println()

Q = bar(distY; label="", ylim=(0, 0.62), c=2)
title!("distribution of Y")

# %%
function show_sim(distX, distY, m, n; L=10^5, P=P, Q=Q)
    ecdf_pval = sim(distX, distY, m, n; L)

    R = plot(ecdf_pval, 0, 0.1; label="")
    plot!(identity; label="", ls=:dot, c=:gray)
    plot!(xguide="α", yguide="probability of p-value ≤ α")
    plot!(xtick=0:0.01:1, ytick=0:0.01:1)
    title!("Mann-Whitney U-test (sizes of X, Y = $m, $n)")
    plot!(size=(400, 400))

    @show distX distY
    @show mean(distX), var(distX), kurtosis(distX), mu(6, distX)
    @show mean(distY), var(distY), kurtosis(distY), mu(6, distY)
    println()
    @show ecdf_pval(0.05)
    println()

    plot(P, Q, R; layout=@layout[[a; b] c], size=(800, 400))
end

# %%
show_sim(distX, distY, 50, 50)

# %%
show_sim(distX, distY, 50, 100)

# %%
show_sim(distX, distY, 50, 200)

# %%
