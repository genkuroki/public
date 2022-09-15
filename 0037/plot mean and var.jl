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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, size=(500, 300),
    titlefontsize=10, tickfontsize=6, guidefontsize=10)

distname(dist) = replace(string(dist), r"{[^}]*}"=>"")

function distname(dist::LocationScale)
    μ, σ, ρ = params(dist)
    if μ == 0
        "$σ * $(distname(ρ))"
    else
        "$μ + $σ * $(distname(ρ))"
    end
end

function plot_mean_and_var(; dist = Normal(10, 2), n = 20, L = 10^6, M = 10^4,
        binX̄=:auto, binS²=:auto, xlimX̄=:auto, xlimS²=:auto, kwargs...)
    X̄ = Vector{Float64}(undef, L)
    S² = Vector{Float64}(undef, L)
    tmpX = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, tmpX[Threads.threadid()])
        X̄[i] = mean(X)
        S²[i] = var(X)
    end
    
    @show μ = mean(dist)
    @show σ² = var(dist)
    @show κ₃ = skewness(dist)
    @show κ₄ = kurtosis(dist)
    @show varS² = σ²^2*(κ₄/n + 2/(n-1))
    @show a = varS²/(2σ²)
    @show ν = 2σ²^2/varS²
    @show normal = Normal(μ, √(σ²/n))
    achisqν = a*Chisq(ν)
    println("achisqν = a * Chisq(ν) = ", distname(achisqν))
    @show normalS² = Normal(σ², √varS²)
    
    P1 = scatter(X̄[1:M], S²[1:M]; label="", ms=1, msw=0, alpha=0.5)
    plot!(xguide="\$\\overline{X}\$", yguide="\$S^2\$", yscale=:log10)
    title!("$(distname(dist)), n=$n")
    plot!(; kwargs...)
    P2 = stephist(X̄; norm=true, bin=binX̄, xlim=xlimX̄, 
        label="", title="\$\\overline{X}\$")
    plot!(normal, quantile.(Ref(X̄), (0.000005, 0.999995))...; label="", ls=:dash)
    P3 = stephist(S²; norm=true, bin=binS², xlim=xlimS², 
        label="", title="\$S^2\$")
    xmax = min(quantile(S², ν ≤ 2 ? 0.999 : 0.99995))
    if binS² != :auto
        xmax = min(last(binS²), xmax)
    end
    if xlimS² != :auto
        xmax = min(last(xlimS²), xmax)
    end
    xmin = max(quantile(S², 0.000005), xmax/200)
    plot!(normalS², xmin, xmax; label="", ls=:dash)
    plot!(achisqν, xmin, xmax; label="", ls=:dashdot)
    
    layout = @layout [a [b; c]]
    plot(P1, P2, P3; layout, size=(800, 400))
    plot!(leftmargin=4Plots.mm)
end

# %%
plot_mean_and_var(dist = Normal(10, 2), n = 20)

# %%
plot_mean_and_var(dist = Uniform(), n = 20)

# %%
plot_mean_and_var(dist = Gamma(4, 2), n = 20, xlimS² = (0, 60))

# %%
plot_mean_and_var(dist = Poisson(1), n = 20, binX̄=0:0.05:2, binS²=0:0.15:3)

# %%
plot_mean_and_var(dist = Poisson(1), n = 100, binX̄=0.6:0.01:1.5, binS²=0.3:0.05:2)

# %%
plot_mean_and_var(dist = LogNormal(), n = 20, xlimX̄=(0.5, 4), binS²=0:0.05:20)

# %%
plot_mean_and_var(dist = LogNormal(), n = 100, xlimX̄=(0.9, 2.4), binS²=0:0.05:20)

# %%
plot_mean_and_var(dist = LogNormal(), n = 500, xlimX̄=(1.2, 2.2), xlimS²=(1, 12))

# %%
@time plot_mean_and_var(dist = LogNormal(), n = 2500, xlimS²=(2.5, 10))

# %%
@time plot_mean_and_var(dist = LogNormal(), n = 12500, xlimS²=(3, 7))

# %%
@time plot_mean_and_var(dist = LogNormal(), n = 62500, xlimS²=(4, 5.7))

# %%
