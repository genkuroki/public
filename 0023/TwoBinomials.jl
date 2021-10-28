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
using Distributions

struct TwoBinomials{T<:Real} <: DiscreteMultivariateDistribution
    m::Int
    p::T
    n::Int
    q::T
end

Base.length(d::TwoBinomials) = 4
Base.eltype(d::TwoBinomials) = Int
Distributions.params(d::TwoBinomials) = (d.m, d.p, d.n, d.q)

function Distributions._rand!(rng::Distributions.AbstractRNG, d::TwoBinomials, x::AbstractVector)
    k = rand(rng, Binomial(d.m, d.p))
    l = rand(rng, Binomial(d.n, d.q))
    @inbounds x[1], x[2], x[3], x[4] = k, l, d.m - k, d.n - l
    x
end

function Distributions._rand!(rng::Distributions.AbstractRNG, d::TwoBinomials, x::AbstractMatrix)
    @inbounds for j in axes(x, 2) 
        k = rand(rng, Binomial(d.m, d.p))
        l = rand(rng, Binomial(d.n, d.q))
        x[1,j], x[2,j], x[3,j], x[4,j] = k, l, d.m - k, d.n - l
    end
    x
end

function Distributions._logpdf(d::TwoBinomials, x::AbstractVector)
    x[1] + x[3] == d.m || x[2] + x[4] == d.m || return -Inf
    logpdf(Binomial(d.m, d.p), x[1]) + logpdf(Binomial(d.n, d.q), x[2]) 
end

Distributions.mean(d::TwoBinomials) = [d.m*d.p, d.n*d.q, d.m*(1 - d.p), d.n*(1 - d.q)]

# %%
d = TwoBinomials(10, 0.25, 20, 0.75)

# %%
params(d)

# %%
rand(d)

# %%
rand(d, 16)

# %%
pdf(d, [3, 16, d.m-3, d.n-16])

# %%
mean(d)

# %%
