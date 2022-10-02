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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://twitter.com/julia_kizi/status/1576613736755912704

# %%
using Distributions
using Random
using StatsPlots

dice = '⚀':'⚅'
ntrials = 10^4
ss = samplesize = 6
tmp = Vector{Float64}(undef, samplesize)
dist = Categorical(0.5, 0.01, 0.01, 0.08, 0.4, 0.0)
means = [mean(rand!(dist, tmp)) for _ in 1:ntrials]
histogram(means; norm=true, bins=10, label="")
title = "$((dist.p...,)), ss=$ss"
plot!(; title)

# %%
using Distributions
using Random
using StatsPlots

dice = '⚀':'⚅'
ntrials = 10^4
ss = samplesize = 6
tmp = Vector{Float64}(undef, samplesize)
dist = Categorical(0.5, 0.01, 0.01, 0.08, 0.4, 0.0)
means = [mean(rand!(dist, tmp)) for _ in 1:ntrials]
histogram(means; norm=true, bin=1-1/(2ss):1/ss:6+1/(2ss), label="")
title = "$((dist.p...,)), ss=$ss"
plot!(; title)

# %%
using Distributions
using Random
using StatsPlots

dice = '⚀':'⚅'
ntrials = 1000
ss = samplesize = 6
dist = Categorical(0.5, 0.01, 0.01, 0.08, 0.4, 0.0)
Xs = [rand(dist, ss) for _ in 1:ntrials]
means = mean.(Xs)
@time anim = @animate for t in [1:ntrials; fill(ntrials, 40)]
    X = Xs[t]
    @views histogram(means[1:t]; norm=true, bin=1-1/(2ss):1/ss:6+1/(2ss), label="")
    title = "$((dist.p...,)), ss=$ss\n"
    title *= "$t times: "
    title *= "($(dice[X[1]])" * prod("+$(dice[X[i]])" for i in 2:ss) * ")/6"
    plot!(; title, ylim=(0.0, 1.2))
end

gif(anim, "result6.gif")

# %%
using Distributions
using Random
using StatsPlots

dice = '⚀':'⚅'
ntrials = 1000
ss = samplesize = 6
dist = Categorical(fill(1//6, ss))
Xs = [rand(dist, ss) for _ in 1:ntrials]
means = mean.(Xs)
@time anim = @animate for t in [1:ntrials; fill(ntrials, 40)]
    X = Xs[t]
    @views histogram(means[1:t]; norm=true, bin=1-1/(2ss):1/ss:6+1/(2ss), label="")
    title = "$((dist.p...,)), ss=$ss\n"
    title *= "$t times: "
    title *= "($(dice[X[1]])" * prod("+$(dice[X[i]])" for i in 2:ss) * ")/6"
    plot!(; title, ylim=(0.0, 1.2))
end

gif(anim, "result.gif")

# %%
