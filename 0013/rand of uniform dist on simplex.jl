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
#     display_name: Julia 1.7.0-beta3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Random, Statistics, Plots

function diff!(X)
    @inbounds for i in reverse(keys(X)[begin:end-1])
        X[i+1] -= X[i]
    end
    X
end

function randsimplex!(rng::AbstractRNG, X)
    X[end] = one(eltype(X))
    Y = @view X[begin:end-1]
    rand!(rng, Y)
    sort!(Y)
    X .*= length(X)
    diff!(X)
end
randsimplex(rng::AbstractRNG, N, T=Float64) = randsimplex!(rng, Vector{T}(undef, N))
randsimplex!(X) = randsimplex!(Random.default_rng(), X)
randsimplex(N, T=Float64) = randsimplex!(Vector{T}(undef, N))

n = 10^4
A = randsimplex(n)
@show mean(A), var(A)
histogram(A; norm=true, alpha=0.3, label="randsimplex($n)")
plot!(x->exp(-x), 0, 6; xlim=(-0.1, 6), lw=2, label="Exponential(1) dist")

# %%
struct Simplex{T} N::Int end
Simplex(N) = Simplex{Float64}(N)

function Random.rand!(rng::AbstractRNG,
        X::AbstractVector{T}, d::Random.SamplerTrivial{Simplex{T}}) where T
    @assert length(X) == d[].N
    randsimplex!(rng, X)
end
Base.rand(rng::AbstractRNG, d::Random.SamplerTrivial{Simplex{T}}) where T =
    rand!(rng, Vector{T}(undef, d[].N), d)

# %%
using BenchmarkTools

X = zeros(10)
d = Simplex(10)
@btime rand!($X, $d)

# %%
rand(Simplex(10))

# %%
using Random, LinearAlgebra, Statistics, Plots

function randsphere!(rng::AbstractRNG, X)
    randn!(X)
    X .*= √(length(X)) / norm(X)
end
randsphere(rng::AbstractRNG, N, T=Float64) = randsphere!(rng, Vector{T}(undef, N))
randsphere!(X) = randsphere!(Random.default_rng(), X)
randsphere(N, T=Float64) = randsphere!(Vector{T}(undef, N))

n = 10^4
A = randsphere(n)
@show mean(A), var(A)
histogram(A; norm=true, alpha=0.3, xlim=(-5, 5), label="randsphere($n)")
plot!(x->exp(-x^2/2)/√(2π), -5, 5; lw=2, label="standard normsl dist")

# %%
