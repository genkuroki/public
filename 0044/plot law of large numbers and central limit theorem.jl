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
using StatsPlots
default(fmt=:png)

dist = Gamma(2, 3)
μ, σ = mean(dist), std(dist)
n = 300

L = 10^6
X̄ = zeros(L)
Threads.@threads for i in 1:L
    X = rand(dist, n)
    X̄[i] = mean(X)
end

P = stephist(X̄; norm=true, label="")
plot!(xlim=(μ-5σ, μ+5σ))
vline!([μ]; label="μ", ls=:dot, c=3)
plot!(ytick=[0])

Q = stephist(X̄; norm=true, label="")
plot!(xlim=(μ-5σ/√n, μ+5σ/√n))
plot!(Normal(μ, σ/√n); label="Normal(μ, σ/√n)", ls=:dash)
vline!([μ]; label="μ", ls=:dot)
plot!(ytick=[0])

plot(P, Q; size=(600, 500), layout=(2, 1))
plot!(plot_title="distribution of sample means for n = $n")

# %%
using Distributions
using StatsPlots
default(fmt=:png)

dist = Poisson(1)
dist = Bernoulli(0.1)
m = 5

μ, σ = mean(dist), std(dist)
n = 1000

L = 10^6
X̄ = zeros(L)
Threads.@threads for i in 1:L
    X = rand(dist, n)
    X̄[i] = mean(X)
end

nbin = 101
binmin, binmax = floor(n*(μ-m*σ/√n))/n, ceil(n*(μ+m*σ/√n))/n
binstep = max(1, round(n*(binmax - binmin)/nbin))/n
bin = binmin-binstep/2:binstep:binmax+binstep/2

P = stephist(X̄; norm=true, bin, label="")
plot!(xlim=(μ-m*σ, μ+m*σ))
vline!([μ]; label="μ", ls=:dot, c=3)
plot!(ytick=[0])

Q = stephist(X̄; norm=true, bin, label="")
plot!(xlim=(μ-m*σ/√n, μ+m*σ/√n))
plot!(Normal(μ, σ/√n); label="Normal(μ, σ/√n)", ls=:dash)
vline!([μ]; label="μ", ls=:dot)
plot!(ytick=[0])

plot(P, Q; size=(600, 500), layout=(2, 1))
plot!(plot_title="distribution of sample means for n = $n")

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

function plot_lln_clt(; dist=Gamma(2, 3), m=5, nbin=200, L=10^6)
    μ, σ = mean(dist), std(dist)
    ns = round.(Int, exp.(range(log(1), log(1000), 200)))
    nmax = maximum(ns)

    X̄ = zeros(L, length(ns))
    Xtmp = [zeros(eltype(dist), nmax) for _ in 1:Threads.nthreads()]
    @time Threads.@threads for i in 1:L
        X = rand!(dist, Xtmp[Threads.threadid()])
        for (j, n) in enumerate(ns)
            X̄[i, j] = mean(@view(X[1:n]))
        end
    end

    js = collect(eachindex(ns))
    js_gif = [fill(js[begin], 60); js; fill(js[end], 60)]
    ns_gif = ns[js_gif]
    @time @gif for (j, n) in zip(js_gif, ns_gif)
        bin = if dist isa DiscreteUnivariateDistribution
            binmin, binmax = floor(n*(μ-m*σ/√n))/n, ceil(n*(μ+m*σ/√n))/n
            binstep = max(1, round(n*(binmax - binmin)/nbin))/n
            binmin-binstep/2:binstep:binmax+binstep/2
        else
            :auto
        end
        
        P = stephist(@view(X̄[:, j]); norm=true, bin, label="")
        plot!(xlim=(μ-m*σ, μ+m*σ))
        vline!([μ]; label="μ", ls=:dot, c=3)
        title!("law of large numbers (n = $n)"; titlefontsize=12)
        plot!(ytick=[0])

        Q = stephist(@view(X̄[:, j]); norm=true, bin, label="")
        plot!(xlim=(μ-m*σ/√n, μ+m*σ/√n))
        plot!(Normal(μ, σ/√n), μ-m*σ/√n, μ+m*σ/√n; label="Normal(μ, σ/√n)", ls=:dash)
        vline!([μ]; label="μ", ls=:dot)
        title!("central limit theorem (n = $n)"; titlefontsize=12)
        plot!(ytick=[0])

        plot(P, Q; size=(600, 500), layout=(2, 1))
        plot!(plot_title="distribution of sample means for n = $n",
            plot_titlefontsize=14)
    end
end

# %%
plot_lln_clt(; dist = Uniform())

# %%
plot_lln_clt(; dist = Gamma(2, 3))

# %%
plot_lln_clt(; dist = MixtureModel([Normal(), Normal(10)], [0.9, 0.1]))

# %%
plot_lln_clt(; dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05]))

# %%
plot_lln_clt(; dist = Poisson(1))

# %%
plot_lln_clt(; dist = Bernoulli(0.1))

# %%
