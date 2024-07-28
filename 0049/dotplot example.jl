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
using StatsPlots
default(fmt=:png)
X = [5652 2187 592 265 13 435 3842 31323 4 500 3842 31323 4 500 352 7 229 284 4 613 883 1556 90 16440 774 2164 776 155 330 10867 4913 2178 16 6488]'
@show mean(X) std(X);
dotplot(X; label="data", size=(320, 400), msc=:auto, ms=3)
hline!([mean(X)]; label="sample mean", legend=:outertop, ls=:dot)
ytick = 0:2000:maximum(X)+2000
plot!(xtick=false, ytick=(ytick, string.(ytick)))

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

X = [5652 2187 592 265 13 435 3842 31323 4 500 3842 31323 4 500 352 7 229 284 4 613 883 1556 90 16440 774 2164 776 155 330 10867 4913 2178 16 6488]' |> collect |> vec

function sim(X, n; L=10^4, samplingfunc! = sample!)
    SampleMean = zeros(L)
    UnbiasedVar = zeros(L)
    tmp = zeros(eltype(X), n)
    for i in 1:L
        Y = samplingfunc!(X, tmp)
        SampleMean[i] = mean(Y)
        UnbiasedVar[i] = var(Y)
    end
    SampleMean, UnbiasedVar
end

function plotsim(X, n; L=10^4, M=10^6, samplingfunc! =sample!)
    μ, σ = mean(X), std(X)
    SampleMean, UnbiasedVar = sim(X, n; L=M, samplingfunc!)
    P = scatter(SampleMean[1:L], UnbiasedVar[1:L]; label="", msc=:auto, ms=1.5, alpha=0.3)
    scatter!([μ], [σ^2]; label="", marker=:star, msc=:auto, c=:red)
    plot!(xguide="sample mean", yguide="sample unbiased variance")
    Q = stephist(SampleMean; norm=true, label="")
    plot!(Normal(μ, σ/√n), extrema(SampleMean)...; label="", ls=:dot)
    R = stephist(UnbiasedVar; permute=(:x, :y), norm=true, label="")
    plot!(σ^2*(Chisq(n-1)/(n-1)), extrema(UnbiasedVar)...; permute=(:x, :y), label="", ls=:dot, yguide=" ")
    plot(Q, P, R; size=(750, 550), layout=@layout[a{0.27h} _; b{0.8w} c], link=:both)
    plot!(plot_title="n = $n")
end

for n in (5, 20, 40, 80, 160, 320)
    plotsim(X, n) |> display
end

# %%
for n in (5, 20, 40, 80, 160, 320)
    plotsim(Normal(2, 3), n; samplingfunc! = rand!) |> display
end

# %%
for n in (5, 20, 40, 80, 160, 320)
    plotsim(Gamma(2, 3), n; samplingfunc! = rand!) |> display
end

# %%
for n in (2, 3, 4, 5, 6, 7, 8)
    plotsim(Uniform(-2, 3), n; samplingfunc! = rand!) |> display
end

# %%
