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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using LinearAlgebra
using Random
using Distributions
using StatsPlots

# %%
A = [
    2 -1
    -1 2
]

# %%
dist = MvNormal(A)

# %%
A[1,2]/√(A[1,1]*A[2,2])

# %%
n = 2^10
XY = rand(dist, n)
X, Y = XY[1,:], XY[2,:]
@show cor(X, Y)

scatter(X, Y)

# %%
function sim(A, n, L=10^5)
    dist = MvNormal(A)
    R = zeros(L)
    XY = zeros(2, n)
    for i in eachindex(R)
        rand!(dist, XY)
        R[i] = @views cor(XY[1,:], XY[2,:])
    end
    R
end

# %%
n = 2^6
R = sim(A, n)
Z = atanh.(R)
r0 = A[1,2]/√(A[1,1]*A[2,2])
z0 = atanh(r0)
histogram(Z; norm=true, alpha=0.3, bin=200)
plot!(Normal(z0, 1/√(n - 3)))
plot!(fit_mle(Normal, Z))

# %%
n = 2^8
x, y = rand(Uniform(-1, 1), n), rand(Uniform(-1, 1), n)
X, Y = x, x + y
@show cor(X, Y)
scatter(X, Y; label="")

# %%
function sim2(n, L=10^5)
    R = zeros(L)
    d = Uniform(-1, 1)
    x = zeros(n)
    y = similar(x)
    Y = similar(x)
    for i in eachindex(R)
        rand!(d, x)
        rand!(d, y)
        @. Y = x + y
        R[i] = cor(x, Y)
    end
    R
end

# %%
n = 2^6
R = sim2(n)
Z = atanh.(R)
r0 = 1/√2
z0 = atanh(r0)
histogram(Z; norm=true, alpha=0.3, bin=200)
plot!(Normal(z0, 1/√(n - 3)))
plot!(fit_mle(Normal, Z))

# %%
n = 2^8
x, y = rand(Exponential(), n), rand(Exponential(), n)
X, Y = x, x + y
@show cor(X, Y)
scatter(X, Y; label="")

# %%
function sim3(n, L=10^5)
    R = zeros(L)
    d = Exponential()
    x = zeros(n)
    y = similar(x)
    Y = similar(x)
    for i in eachindex(R)
        rand!(d, x)
        rand!(d, y)
        @. Y = x + y
        R[i] = cor(x, Y)
    end
    R
end

# %%
n = 2^6
R = sim3(n)
Z = atanh.(R)
r0 = 1/√2
z0 = atanh(r0)
histogram(Z; norm=true, alpha=0.3, bin=200)
plot!(Normal(z0, 1/√(n - 3)))
plot!(fit_mle(Normal, Z))

# %%
function simplot(dist, n, L=10^5)
    Z = zeros(L)
    x = zeros(n)
    y = similar(x)
    Y = similar(x)
    for i in eachindex(R)
        rand!(dist, x)
        rand!(dist, y)
        @. Y = x + y
        Z[i] = atanh(cor(x, Y))
    end
    r0 = 1/√2
    z0 = atanh(r0)
    histogram(Z; norm=true, alpha=0.3, bin=200)
    plot!(Normal(z0, 1/√(n - 3)))
    plot!(fit_mle(Normal, Z))
end

# %%
simplot(Normal(), 2^6)

# %%
simplot(Uniform(), 2^6)

# %%
simplot(Exponential(), 2^6)

# %%
simplot(MixtureModel([Normal(), Normal(10, 1)], [0.95, 0.05]), 2^4)

# %%
simplot(MixtureModel([Normal(), Normal(10, 1)], [0.95, 0.05]), 2^6)

# %%
