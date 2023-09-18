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
using QuadGK
using StatsPlots
default(fmt=:png)

# %%
Normal()

# %%
2 + 3Normal()

# %%
plot(2 + 3Normal())

# %%
Gamma(2, 3)

# %%
100 + 5Gamma(2, 3)

# %%
plot(100 + 5Gamma(2, 3))

# %%
plot(10 + 2TDist(4); label="10 + 2TDist(4)")
plot!(Normal(10, 4); label="Normal(10, 4)")

# %%
prod_poissons = product_distribution([Poisson(7), Poisson(3)])

# %%
X = rand(prod_poissons, 10^6)

# %%
Y = [v[1] for v in eachcol(X) if sum(v) == 10]
stephist(Y; norm=true, bin=-0.5:10.5, xtick=0:10,
    label="Poisson(7)×Poisson(3) | sum = 10")
bar!(Binomial(10, 7/(7+3)); alpha=0.3, label="Binomial(10, 7/(7+3))")

# %%
mixnormal = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])

# %%
plot(mixnormal) # maybe bug

# %%
plot(1mixnormal)

# %%
1mixnormal

# %%
L = 10^6
X = rand(Gamma(7, 2), L)
Y = rand(Gamma(3, 2), L)
Z = @. X / (X + Y)

stephist(Z; norm=true, label="X/(X+Y) where X~Gamma(7,2), Y~Gamma(3,2)")
plot!(Beta(7, 3); label="Beta(7, 3)", ls=:dash)
plot!(legend=:outertop)

# %%
L = 10^6
X = rand(Gamma(7, 2), L)
Y = rand(Gamma(3, 2), L)
Z = @. X / Y

stephist(Z; norm=true, label="X/Y where X~Gamma(7,2), Y~Gamma(3,2)")
plot!(BetaPrime(7, 3); label="BetaPrime(7, 3)", ls=:dash)
plot!(legend=:outertop, xlim=(-0.1, 20))

# %%
L = 10^6
X = rand(Binomial(10, 0.5), L)
Y = rand(Binomial(10, 0.5), L)
Z = X + Y

@show var(Z)
@show var(Binomial(20, 0.5))

stephist(Z; norm=true, bin=-0.5:20.5, label="Binomial(10, 0.5) + Binomial(10, 0.5)")
bar!(Binomial(20, 0.5); alpha=0.3, label="Binomial(20, 0.5)")
plot!(legend=:outertop)

# %%
L = 10^6
X = rand(Binomial(10, 0.2), L)
Y = rand(Binomial(10, 0.8), L)
Z = X + Y

@show var(Z)
@show var(Binomial(20, 0.5))
@show var(Binomial(10, 0.2)) + var(Binomial(10, 0.8))

stephist(Z; norm=true, bin=-0.5:20.5, label="Binomial(10, 0.2) + Binomial(10, 0.8)")
bar!(Binomial(20, 0.5); alpha=0.3, label="Binomial(20, 0.5)")
plot!(Normal(10, √3.2); label="Normal(10, √3.2)", ls=:dash)
plot!(legend=:outertop, ylim=(-0.005, 0.25))

# %%
