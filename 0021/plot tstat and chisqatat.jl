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
#     display_name: Julia 1.6.3
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
    plot!(TDist(n-1); lw=1.5, label="TDist(n - 1)")
    sk = skewness(dist)
    diststr = replace("$dist", r"\{.*\}"=>"")
    title!("$diststr,  n = $n\nskewness = $(round(sk; digits=3))"; titlefontsize=10)
end

function plot_chisqstat(;
        n = 10,
        L = 10^5,
        dist = Normal(1, 2),
    )
    X = (n-1)*[var(rand(dist, n)) for _ in 1:L]/var(dist)
    m, s = mean(X), std(X)
    xlim = (max(-1, m - 4s), n + 4s)
    bin = range(0, xlim[2]; length=round(Int, 51/(1 - xlim[1]/xlim[2])))
    histogram(X; norm=true, alpha=0.3, label="(n - 1)u²/σ²", xlim, bin)
    plot!(Chisq(n-1); lw=1.5, label="Chisq(n - 1)")
    kurt = kurtosis(dist)
    if isfinite(kurt)
        M, S = n - 1, (n - 1)*√((kurt + 3)/n - (n - 3)/(n*(n - 1)))
        plot!(Normal(M, S); lw=1.5, ls=:dash, label="normal dist.")
    end
    diststr = replace("$dist", r"\{.*\}"=>"")
    title!("$diststr,  n = $n\nkurtosis = $(round(kurt; digits=3))"; titlefontsize=10)
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
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; n) |> display
end

# %%
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; dist = Uniform(0, 1), n) |> display
end

# %%
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; dist = TDist(4), n) |> display
end

# %%
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; dist = TDist(5), n) |> display
end

# %%
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; dist = Gamma(2, 1), n) |> display
end

# %%
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; dist = Gamma(4, 1), n) |> display
end

# %%
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; dist = Gamma(10, 1), n) |> display
end

# %%
for n in (10, 20, 50, 100, 200, 400)
    plot_both(; dist = Exponential(), n) |> display
end

# %%
