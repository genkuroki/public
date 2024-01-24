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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using LinearAlgebra
using StatsPlots
default(fmt=:png)

function rand_unit_sphere(n)
    x = randn(n)
    x ./= norm(x)
    x
end

function plot5x6_rand_unit_sphere(; n = 100)
    X = [rand_unit_sphere(n) for _ in 1:30]
    PP = []
    for x in X
        P = stephist(x; norm=true, label="")
        plot!(Normal(0, 1/√n), label="")
        push!(PP, P)
    end
    plot(PP...; size=(1000, 900), layout=(6, 5))
    plot!(plot_title="n = $n")
    plot!(tickfontsize=5)
end

plot5x6_rand_unit_sphere(; n = 100)

# %%
plot5x6_rand_unit_sphere(; n = 1000)

# %%
using Distributions
using LinearAlgebra
using Random
using StatsPlots
default(fmt=:png)

function rand_unit_sphere(n)
    x = randn(n)
    x ./= norm(x)
    x
end

function sim_sum_of_ai_Xi(dist, a; L=10^5)
    Y = zeros(L)
    Xtmp = [similar(a, eltype(dist)) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, Xtmp[Threads.threadid()])
        Y[i] = dot(a, X)
    end
    Y
end

function plot_sum_of_ai_Xi(dist, a; L=10^5, bin=:auto)
    distname = replace(string(dist), r"{[^}]*}"=>"")
    n = length(a)
    μ = sum(a) * mean(dist)
    σ = norm(a) * std(dist)
    @show skewness(dist)
    @show kurtosis(dist)
    @show sum(a.^3) / norm(a)^3
    @show sum(a.^4) / norm(a)^4
    Y = sim_sum_of_ai_Xi(dist, a; L)
    
    P = stephist(collect(a); norm=true, bin, label="", c=3)
    title!("histogram of a with length n=$n")
    
    Q = stephist(Y; norm=true, label="")
    plot!(Normal(μ, σ), label="")
    title!("dist=$distname, sample size n=$n")
    
    plot(P, Q; size=(600, 500), layout=@layout [a{0.3h}; b])
    plot!(titlefontsize=10)
end

n = 100
plot_sum_of_ai_Xi(Exponential(), abs.(randn(n)))

# %%
n = 1000
plot_sum_of_ai_Xi(Exponential(), abs.(randn(n)))

# %%
n = 100
plot_sum_of_ai_Xi(Exponential(), fill(1/n, n); bin=(-0.01:0.02:2.01)/n)

# %%
n = 1000
plot_sum_of_ai_Xi(Exponential(), fill(1/n, n); bin=(-0.01:0.02:2.01)/n)

# %%
n = 100
plot_sum_of_ai_Xi(Exponential(), rand(n))

# %%
n = 1000
plot_sum_of_ai_Xi(Exponential(), rand(n))

# %%
n = 100
plot_sum_of_ai_Xi(Exponential(), rand(n).^2)

# %%
n = 1000
plot_sum_of_ai_Xi(Exponential(), rand(n).^2)

# %%
n = 100
plot_sum_of_ai_Xi(Exponential(), rand(n).^3)

# %%
n = 1000
plot_sum_of_ai_Xi(Exponential(), rand(n).^3)

# %%
