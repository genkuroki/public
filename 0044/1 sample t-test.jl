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
@time using StatsPlots
default(fmt=:png)

# %%
X = [fill(1, 25); fill(2, 12); fill(3, 3); fill(4, 3); fill(5, 3); 22]
X = X - rand(length(X))
@show X;

# %%
histogram(X)

# %%
(1 - 1/length(X))^length(X), 1/exp(1)

# %%
function _sim(f!, μ, n, t; L=10^6)
    T = Vector{Float64}(undef, L)
    Ytmp = [Vector{t}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        Y = f!(Ytmp[Threads.threadid()])
        t = (mean(Y) - μ)/√(var(Y)/n)
        T[i] = t
    end
    T
end

# %%
function sim(X::AbstractVector; L=10^6)
    μ = mean(X)
    n = length(X)
    f!(Ytmp) = sample!(X, Ytmp)
    _sim(f!, μ, n, eltype(X); L)
end

# %%
function sim(dist::UnivariateDistribution, n; L=10^6)
    μ = mean(dist)
    f!(Ytmp) = rand!(dist, Ytmp)
    _sim(f!, μ, n, eltype(dist); L)
end

# %%
T = sim(X)
stephist(T; norm=true, bin=2^7)
plot!(Normal())

# %%
T = sim(X)
stephist(T .^ 2; norm=true)
plot!(Chisq(1))
plot!(xlim=(-0.2, 6), ylim=(-0.02, 1))

# %%
T = sim(X[1:end-1])
stephist(T; norm=true, bin=2^7)
plot!(Normal())

# %%
T = sim(X[1:end-1])
stephist(T .^ 2; norm=true)
plot!(Chisq(1))
plot!(xlim=(-0.2, 6), ylim=(-0.02, 1))

# %%
T = sim(Exponential(), 31)
stephist(T; norm=true, bin=2^7)
plot!(Normal())

# %%
T = sim(Exponential(), 31)
stephist(T .^ 2; norm=true)
plot!(Chisq(1))
plot!(xlim=(-0.2, 6), ylim=(-0.02, 1))

# %%
