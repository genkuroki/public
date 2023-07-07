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
        X = rand!(dist, Xtmp[tid])
        X̄ = mean(X)
        S = std(X)
        T = √n * (X̄ - μ) / S
        pval[i] = 2ccdf(D, abs(T))
    end
    _ecdf = ecdf(pval)
    f(x) = _ecdf(x)
    println("n = $n => coverage probability (α=5%) = ", 1 - f(0.05))
    plot(f, 0, 0.1; label="")
    plot!(identity; label="", c=:grey, ls=:dash)
    plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
    plot!(xguide="α", yguide="probability of P-value ≤ α")
    title!("$testname$distname, n=$n", titlefontsize=10)
    plot!(size=(400, 420))
end

function plot_sims(dist; ns=(5, 10, 20, 100), correction=true)
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
