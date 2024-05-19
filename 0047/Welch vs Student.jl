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
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using HypothesisTests
using Random
using StatsPlots
default(fmt=:png)

# %%
bar(Poisson(5); label="Poisson(5)", alpha=0.3)
plot!(Normal(5, √5); label="Normal(5, √5)")

# %%
bar(Poisson(0.05); label="Poisson(0.05)", alpha=0.3)
plot!(Normal(0.05, √0.05); label="Normal(0.05, √0.05)")
plot!(xtick=-10:100)

# %%
bar(Poisson(0.10); label="Poisson(0.10)", alpha=0.3)
plot!(Normal(0.10, √0.10); label="Normal(0.10, √0.10)")
plot!(xtick=-10:100)

# %%
myecdf(A, x) = count(≤(x), A)/length(A)

distname(dist) = replace(string(dist), r"{[^}]*}"=>"")

function distname(dist::LocationScale)
    μ, σ, ρ = params(dist)
    a = μ == 0 ? "" : "$μ + "
    b = σ == 1 ? "" : "$σ "
    a * b * distname(ρ)
end

function plot_ecdf_pvals(; distx=Poisson(0.05), disty=Poisson(0.05), m=100, n=900, L=10^5)
    pval_welch = zeros(L)
    pval_student = zeros(L)
    nth = Threads.nthreads()
    Xtmp = [zeros(m) for _ in 1:nth]
    Ytmp = [zeros(n) for _ in 1:nth]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(distx, Xtmp[tid])
        Y = rand!(disty, Ytmp[tid])
        if mean(X) == mean(Y)
            pval_welch[i] = 1.0
            pval_student[i] = 1.0
        else
            pval_welch[i] = pvalue(UnequalVarianceTTest(X, Y))
            pval_student[i] = pvalue(EqualVarianceTTest(X, Y))
        end
    end
    
    println("distx = ", distname(distx))
    println("disty = ", distname(disty))
    @show m n
    @show mean(distx) - mean(disty)
    @show std(distx) / std(disty)

    plot(α->myecdf(pval_welch, α), 0, 0.1; label="Welch")
    plot!(α->myecdf(pval_student, α); label="Student")
    plot!(identity; label="", ls=:dot, c=:gray)
    plot!(xtick=0:0.01:1, ytick=0:0.01:1)
    plot!(xguide="α", yguide="probability of P-value ≤ α")
    plot!(size=(400, 400))
end

# %%
plot_ecdf_pvals(; distx=Poisson(0.05), disty=Poisson(0.05), m=100, n=900)

# %%
plot_ecdf_pvals(; distx=Poisson(0.10), disty=Poisson(0.05)+0.05, m=100, n=900)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05)+0.05, disty=Poisson(0.10), m=100, n=900)

# %%
plot_ecdf_pvals(; distx=Poisson(1.2^2*0.05), disty=Poisson(0.05)+(1.2^2-1)*0.05, m=100, n=900)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05)+(1.2^2-1)*0.05, disty=Poisson(1.2^2*0.05), m=100, n=900)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05), disty=Poisson(0.05), m=100, n=100)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05), disty=Poisson(0.05), m=100, n=200)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05), disty=Poisson(0.05), m=200, n=1800)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05), disty=Poisson(0.05), m=400, n=3600)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05), disty=Poisson(0.05), m=1000, n=9000)

# %%
plot_ecdf_pvals(; distx=Poisson(0.10), disty=Poisson(0.05)+0.05, m=1000, n=9000)

# %%
plot_ecdf_pvals(; distx=Poisson(0.05)+0.05, disty=Poisson(0.10), m=1000, n=9000)

# %%

# %%
plot_ecdf_pvals(; distx=Normal(0, 1), disty=Normal(0, 4), m=200, n=100)

# %%
plot_ecdf_pvals(; distx=Normal(0, 1), disty=Normal(0, 4), m=100, n=200)

# %%

# %%
dist = InverseGamma(3.1, 5)
dist = dist - mean(dist)
plot(dist; label="dist")
plot!(xlim=(-3, 8))

# %%
plot_ecdf_pvals(; distx=dist, disty=dist, m=80, n=20)

# %%
plot_ecdf_pvals(; distx=dist, disty=dist, m=320, n=80)

# %%
plot_ecdf_pvals(; distx=dist, disty=2dist, m=80, n=20)

# %%
plot_ecdf_pvals(; distx=dist, disty=2dist, m=20, n=80)

# %%
plot_ecdf_pvals(; distx=dist, disty=2dist, m=320, n=80)

# %%
plot_ecdf_pvals(; distx=dist, disty=2dist, m=80, n=320)

# %%
