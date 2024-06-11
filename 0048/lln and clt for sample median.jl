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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

function plot_lln_clt_median(dist, n; L=10^5)
    @show dist
    @show n
    @show m = median(dist) # population median
    @show s = 1/(2pdf(dist, m))
    @show se = s/√n # std of sample median
    @show normal = Normal(m, se)
    M = zeros(L)
    Xtmp = [zeros(n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid])
        M[i] = median(X)
    end
    P = stephist(M; norm=true, label="sample medians")
    plot!(xlim=(m - 5s, m + 5s))
    vline!([m]; label="population median", ls=:dash)
    title!("law of large numbers for median")
    Q = stephist(M; norm=true, label="sample medians")
    plot!(normal; label="normal approx.", ls=:dash)
    plot!(xlim=(m - 5se, m + 5se))
    title!("central limit theorem for median")
    plot(P, Q; size=(600, 600), layout=(2, 1))
end

# %%
plot_lln_clt_median(Exponential(), 10)

# %%
plot_lln_clt_median(Exponential(), 30)

# %%
plot_lln_clt_median(Exponential(), 100)

# %%
plot_lln_clt_median(Exponential(), 300)

# %%
dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
plot(dist, -25, 25; label="dist", size=(600, 300))

# %%
plot_lln_clt_median(dist, 10)

# %%
plot_lln_clt_median(dist, 30)

# %%
plot_lln_clt_median(dist, 100)

# %%
plot_lln_clt_median(dist, 300)

# %%
