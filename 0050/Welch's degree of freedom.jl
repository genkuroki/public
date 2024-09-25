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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, titlefontsize=12)

function degree_of_freedom_welch(m, sx², n, sy²)
    (sx²/m + sy²/n)^2 / ((sx²/m)^2/(m-1) + (sy²/n)^2/(n-1))
end

function degree_of_freedom_welch(x, y)
    m, sx² = length(x), var(x)
    n, sy² = length(y), var(y)
    degree_of_freedom_welch(m, sx², n, sy²)
end

function plot_df_welch(; distx=Normal(0, 1), disty=Normal(0, 1), m=5, n=5, L=10^7)
    xtmp, ytmp = zeros(m), zeros(n)
    df_welch = [degree_of_freedom_welch(rand!(distx, xtmp), rand!(disty, ytmp)) for _ in 1:L]
    println("m = ", m, ", n = ", n, ", extrema(df_welch) = ", extrema(df_welch))
    stephist(df_welch; norm=true, label="df_welch")
    title!("var(distx)=$(var(distx)), var(disty)=$(var(disty)), m=$m, n=$n")
end

PP = []
for (m, n) in ((5, 5), (10, 10), (5, 10), (5, 15))
    P = plot_df_welch(; distx=Normal(0, 1), disty=Normal(0, 1), m, n)
    push!(PP, P)
end
println()
plot(PP...; size=(1000, 600), layout=(2, 2))

# %%
PP = []
for (m, n) in ((5, 5), (10, 10), (5, 10), (5, 15), (10, 5), (15, 5))
    P = plot_df_welch(; distx=Normal(0, 1), disty=Normal(0, 2), m, n)
    push!(PP, P)
end
println()
plot(PP...; size=(1000, 900), layout=(3, 2))

# %%
