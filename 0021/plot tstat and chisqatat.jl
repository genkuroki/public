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

# %%
using Distributions
using StatsPlots

function plot_tstat(;
        n = 10,
        L = 10^5,
        dist = Normal(1, 2),
        xlim = (-5, 5),
        xtick = -5:5,
        bin = range(extrema(xtick)...; length=51),
    )
    tstat(sample) = (mean(sample) - mean(dist))/√(var(sample)/n)
    T = [tstat(rand(dist, n)) for _ in 1:L]
    histogram(T; norm=true, alpha=0.3, label="(x - μ)/√(u²/n)", xlim, xtick, bin)
    plot!(TDist(n-1); lw=1.5, label="Tdist(n - 1)")
    diststr = replace("$dist", r"\{.*\}"=>"")
    title!("$diststr,  n = $n"; titlefontsize=10)
end

function plot_chisqstat(;
        n = 10,
        L = 10^5,
        dist = Normal(1, 2),
    )
    X = (n-1)*[var(rand(dist, n)) for _ in 1:L]/var(dist)
    m, s = mean(X), std(X)
    xlim = (max(-1, m - 4s), n + 4s)
    bin = range(0, maximum(xlim); length=51)
    histogram(X; norm=true, alpha=0.3, label="(n - 1)u²/σ²", xlim, bin)
    plot!(Chisq(n-1); lw=1.5, label="Chisq(n - 1)")
    diststr = replace("$dist", r"\{.*\}"=>"")
    title!("$diststr,  n = $n"; titlefontsize=10)
end

function plot_both(;
        n = 10,
        L = 10^5,
        dist = Normal(1, 2),
        xlim1 = (-5, 5),
        xtick1 = -5:5,
        bin1 = range(extrema(xtick1)...; length=51),
    )
    P = plot_tstat(; n, L, dist, xlim=xlim1, xtick=xtick1, bin=bin1)
    Q = plot_chisqstat(; n, L, dist)
    plot(P, Q; size=(800, 300))
end

# %%
plot_both()

# %%
plot_both(dist = Uniform(0, 1))

# %%
plot_both(dist = TDist(4))

# %%
plot_both(dist = Gamma(2, 1))

# %%
plot_both(dist = Gamma(4, 1))

# %%
plot_both(dist = Gamma(10, 1))

# %%
plot_both(dist = Exponential())

# %%
plot_both(dist = Exponential(), n = 20)

# %%
plot_both(dist = Exponential(), n = 50)

# %%
plot_both(dist = Exponential(), n = 100)

# %%
plot_both(dist = Exponential(), n = 200)

# %%
plot_both(dist = Exponential(), n = 400)

# %%
