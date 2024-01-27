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
using SpecialFunctions
using StatsPlots
default(fmt=:png)

# %%
module O

using Distributions
using SpecialFunctions: erf, erfinv

struct MyNormal{T<:Real} <: ContinuousUnivariateDistribution
    μ::T
    σ::T
end

MyNormal(μ::Real, σ::Real) = MyNormal(promote(μ, σ)...)
MyNormal(μ::Integer, σ::Integer) = Normal(float(μ), float(σ))

Distributions.params(d::MyNormal) = (d.μ, d.σ)
Distributions.mean(d::MyNormal) = d.μ
Distributions.std(d::MyNormal) = d.σ
Distributions.var(d::MyNormal) = d.σ^2

function Distributions.logpdf(d::MyNormal, x::Real)
    μ, σ = params(d)
    -((x - μ)^2/σ^2 + log(2π*σ^2))/2
end

function Distributions.cdf(d::MyNormal, x::Real)
    μ, σ = params(d)
    (erf((x - μ)/(√2*σ)) + 1)/2
end

function Distributions.quantile(d::MyNormal, p::Real)
    μ, σ = params(d)
    μ + √2 * σ * erfinv(2p - 1)
end

Distributions.maximum(d::MyNormal{T}) where T<:Real = T(Inf)
Distributions.minimum(d::MyNormal{T}) where T<:Real = T(-Inf)

end

# %%
O.MyNormal(3, 2)

# %%
mean(O.MyNormal(3, 2)), std(O.MyNormal(3, 2)), var(O.MyNormal(3, 2))

# %%
logpdf(O.MyNormal(3, 2), 1), logpdf(Normal(3, 2), 1)

# %%
pdf(O.MyNormal(3, 2), 1), pdf(Normal(3, 2), 1)

# %%
logcdf(O.MyNormal(3, 2), 1), logcdf(Normal(3, 2), 1)

# %%
cdf(O.MyNormal(3, 2), 1), cdf(Normal(3, 2), 1)

# %%
quantile(O.MyNormal(3, 2), 0.158655), quantile(Normal(3, 2), 0.158655)

# %%
n = 10^6
X = rand(O.MyNormal(3, 2), n)
#Y = rand(Normal(3, 2), n)

stephist(X; norm=true, label="rand(MyNormal(3, 2), n)")
#stephist!(Y; norm=true, label="rand(Normal(3, 2), n)", ls=:dash)
plot!(Normal(3, 2); label="Normal(3, 2)", ls=:dashdot)

# %%
extrema(O.MyNormal(3, 2))

# %%
insupport(O.MyNormal(3, 2), 1)

# %%
