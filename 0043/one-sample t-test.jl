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
#     display_name: Julia 1.9.2
#     language: julia
#     name: julia-1.9
# ---

# %% tags=[]
using Distributions
using Random
using StatsBase
using StatsPlots

function plot_sim(dist, n; μ=mean(dist), L=10^6, D=TDist(n-1))
    testname = D == TDist(n-1) ? "t-test, " : D == Normal() ? "no correction, " : ""
    distname = replace(string(dist), r"{[^}]*}"=>"")
    pval = Vector{Float64}(undef, L)
    Xtmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = Xtmp[tid]
        rand!(dist, X)
        X̄ = mean(X)
        S = std(X)
        T = √n * (X̄ - μ) / S
        pval[i] = 2ccdf(D, abs(T))
    end
    _ecdf_pval = ecdf(pval)
    ecdf_pval(x) = _ecdf_pval(x)
    println("n = $n => coverage probability (α=5%) = ", 1 - ecdf_pval(0.05))
    plot(ecdf_pval, 0, 0.1; label="")
    plot!(identity; label="", c=:grey, ls=:dash)
    plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
    plot!(xguide="nominl significance level α", yguide="probability of P-value ≤ α")
    title!("$testname$distname, n=$n", titlefontsize=9)
    plot!(size=(400, 420), tickfontsize=6, guidefontsize=8)
end

function plot_sims(dist; ns=(5, 10, 20, 100), correction=true)
    distname = replace(string(dist), r"{[^}]*}"=>"")
    println("skewness of $distname = ", skewness(dist))
    println("kurtosis of $distname = ", kurtosis(dist))
    println()
    PP = []
    for n in ns
        D = correction ? TDist(n-1) : Normal()
        P = plot_sim(dist, n; D)
        push!(PP, P)
    end
    println()
    plot(PP...; size=(800, 840))
end

# %%
plot_sims(Normal(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Normal(); ns=(5, 10, 20, 100))

# %%
plot_sims(Uniform(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Uniform(); ns=(5, 10, 20, 100))

# %%
plot_sims(Exponential(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Exponential(); ns=(5, 10, 20, 100))

# %%
plot_sims(Exponential(); ns=(100, 200, 500, 1000))

# %%
plot_sims(Beta(0.1, 0.1); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(Beta(0.1, 0.1); ns=(5, 10, 20, 100))

# %%
plot_sims(TDist(4); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(TDist(4); ns=(5, 10, 20, 100))

# %%
plot_sims(LogNormal(); ns=(5, 10, 20, 100), correction=false)

# %%
plot_sims(LogNormal(); ns=(5, 10, 20, 100))

# %%
plot_sims(LogNormal(); ns=(100, 300, 1000, 3000), correction=false)

# %%
plot_sims(LogNormal(); ns=(100, 300, 1000, 3000))

# %%

# %%
