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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# https://juliastats.org/Distributions.jl/stable/extends/#Create-a-Distribution

# %%
module O

using Distributions
using SpecialFunctions

struct MyBeta{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
end

Distributions.params(d::MyBeta) = (d.a, d.b)

function Distributions.quantile(d::MyBeta, p::Real)
    a, b = params(d)
    beta_inc_inv(a, b, p)[1]
end

end

@show O.quantile(O.MyBeta(3, 7), 0.5)
@show rand(O.MyBeta(3, 7), 3);

# %%
using Distributions
using StatsPlots

d = O.MyBeta(3, 7)
X = rand(d, 10^5)
histogram(X; norm=true, alpha=0.3, bin=0:0.02:1, label="MyBeta(3, 7)")
plot!(Beta(3, 7); label="Beta(3, 7)")

# %%
