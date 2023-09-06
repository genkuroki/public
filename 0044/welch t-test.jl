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
using Random
using StatsBase
using StatsPlots
default(fmt=:png)

function MyLogNormal(m, s)
    σ² = log(1 + s^2/m^2)
    μ = log(m) - σ²/2
    LogNormal(μ, √σ²)
end

distname(dist) = replace(string(dist), r"{[^\}]*}"=>"")

function distname(dist::LocationScale)
    μ, σ, ρ = params(dist)
    m = μ != 0 ? "$μ + " : ""
    s = σ != 1 ? "$σ " : ""
    m * s * distname(ρ)
end

# %%
function plot_tstats(; distx, m, disty, n, L=10^6, kwargs...)
    Xtmp = [Vector{eltype(distx)}(undef, m) for _ in 1:Threads.nthreads()]
    Ytmp = [Vector{eltype(disty)}(undef, n) for _ in 1:Threads.nthreads()]
    T = Vector{Float64}(undef, L)
    DF = similar(T)
    Pval = similar(T)
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(distx, Xtmp[tid])
        Y = rand!(disty, Ytmp[tid])
        X̄ = mean(X)
        Ȳ = mean(Y)
        sx2 = var(X)
        sy2 = var(Y)
        t = (X̄ - Ȳ)/√(sx2/m + sy2/n)
        T[i] = t
        df = (sx2/m + sy2/n)^2 / ((sx2/m)^2/(m-1) + (sy2/n)^2/(n-1))
        DF[i] = df
        pval = 2ccdf(TDist(df), abs(t))
        Pval[i] = pval
    end
    
    DFmin, DFmax = minimum(DF), maximum(DF)
    _ecdfPval = ecdf(Pval)
    ecdfPval(x) = _ecdfPval(x)
    
    println(distname(distx))
    @show m mean(distx) std(distx)
    println()
    println(distname(disty))
    @show n mean(disty) std(disty)
    println()
    @show DFmin DFmax
    @show ecdfPval(0.05)
    
    P = stephist(T; norm=true, label="t-statistics", lw=1.5, kwargs...)
    plot!(TDist(DFmax); label="", ls=:dash)
    plot!(TDist(DFmin); label="", ls=:dashdot)
    
    T2 = T .^ 2
    _ecdfT2 = ecdf(T2)
    ecdfT2(x) = _ecdfT2(x)
    Q = plot(ecdfT2, 0, 20; label="t²", lw=1.5)
    plot!(x -> cdf(FDist(1, DFmax), x), 0, 20; label="", ls=:dash)
    plot!(x -> cdf(FDist(1, DFmin), x), 0, 20; label="", ls=:dashdot)
    plot!(xlim=(3.5, 8), ylim=(0.9, 1.01))
    plot!(ytick=0:0.01:1)
    
    plot(P, Q; size=(1000, 300))
end

# %%
plot_tstats(distx=Normal(), m=20, disty=Normal(0, 2), n=40)

# %%
plot_tstats(distx=Normal(), m=60, disty=Normal(0, 2), n=20)

# %%
plot_tstats(distx=Normal(), m=40, disty=Normal(0, 2), n=20)

# %%
plot_tstats(distx=LogNormal(), m=30, disty=LogNormal(), n=30)

# %%
plot_tstats(distx=LogNormal(), m=40, disty=LogNormal(), n=20)

# %%
plot_tstats(distx=LogNormal(), m=60, disty=LogNormal(), n=20)

# %%
plot_tstats(distx=LogNormal()-mean(LogNormal()), m=20, disty=2(LogNormal()-mean(LogNormal())), n=40)

# %%
plot_tstats(distx=LogNormal()-mean(LogNormal()), m=30, disty=2(LogNormal()-mean(LogNormal())), n=30)

# %%
plot_tstats(distx=LogNormal()-mean(LogNormal()), m=40, disty=2(LogNormal()-mean(LogNormal())), n=20)

# %%
plot_tstats(distx=LogNormal()-mean(LogNormal()), m=300, disty=2(LogNormal()-mean(LogNormal())), n=300)

# %%
plot_tstats(distx=LogNormal()-mean(LogNormal()), m=400, disty=2(LogNormal()-mean(LogNormal())), n=200)

# %%
plot_tstats(distx=Normal(2), m=20, disty=Normal(2.8), n=40)

# %%
plot_tstats(distx=LogNormal(), m=20, disty=LogNormal()+0.8*std(LogNormal()), n=40)

# %%
plot(100Beta(0.5, 0.5); label="distx")
plot!(100Beta(0.1, 0.1); label="disty", ls=:dash)
plot!(ylim=(-0.003, 0.1))

# %%
plot_tstats(distx=100Beta(0.5, 0.5), m=10, disty=100Beta(0.1, 0.1), n=10)

# %%
plot_tstats(distx=100Beta(0.5, 0.5), m=20, disty=100Beta(0.1, 0.1), n=10)

# %%
distx = MixtureModel([Normal(-2), Normal(2)], [0.5, 0.5])
disty = MixtureModel([Normal(-4), Normal(4)], [0.5, 0.5])

plot(x -> pdf(distx, x), -6, 6; label="distx")
plot!(x -> pdf(disty, x), -10, 10; label="disty", ls=:dash)

# %%
plot_tstats(; distx, m=10, disty, n=10)

# %%
plot_tstats(; distx, m=20, disty, n=10)

# %%
distx = MixtureModel([Normal(-2+1.2), Normal(2+1.2)], [0.8, 0.2])
disty = MixtureModel([Normal(-4+2.4, 2), Normal(4+2.4, 2)], [0.8, 0.2])

plot(x -> pdf(distx, x), -6+1.2, 6+1.2; label="distx")
plot!(x -> pdf(disty, x), -14+2.4, 14+2.4; label="disty", ls=:dash)

# %%
plot_tstats(; distx, m=10, disty, n=10)

# %%
plot_tstats(; distx, m=20, disty, n=20)

# %%
plot_tstats(; distx, m=30, disty, n=30)

# %%
plot_tstats(; distx, m=40, disty, n=40)

# %%
plot_tstats(; distx, m=100, disty, n=100)

# %%
