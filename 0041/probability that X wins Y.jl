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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10)

# %%
probability_that_x_wins_y(x, y) = (x > y) + (x == y)/2

function probability_that_distX_wins_distY(
        distX::DiscreteUnivariateDistribution,
        distY::DiscreteUnivariateDistribution)
    sum(probability_that_x_wins_y(x, y) * pdf(distX, x) * pdf(distY, y)
        for x in support(distX), y in support(distY))
end

# %%
distA = Categorical(0, 1/3, 0, 1/3, 0, 0, 0, 0, 1/3)
distB = Categorical(1/3, 0, 0, 0, 0, 1/3, 0, 1/3, 0)
distC = Categorical(0, 0, 1/3, 0, 1/3, 0, 1/3, 0, 0)

@show probability_that_distX_wins_distY(distA, distB)
@show probability_that_distX_wins_distY(distB, distC)
@show probability_that_distX_wins_distY(distC, distA)
println()

PA = bar(distA; label="", title="distribution of A", c=1)
PB = bar(distB; label="", title="distribution of B", c=2)
PC = bar(distC; label="", title="distribution of C", c=3)
plot(PA, PB, PC; layout=(3, 1))

# %%
distX = Categorical(0, 2/3, 0, 0, 1/3)
distY = Categorical(1/3, 0, 0, 2/3, 0)
distZ = Categorical(0, 1/3, 1/3, 1/3, 0)

@show probability_that_distX_wins_distY(distX, distY)
@show probability_that_distX_wins_distY(distY, distZ)
@show probability_that_distX_wins_distY(distZ, distX)
println()

PX = bar(distX; label="", title="distribution of X", c=1)
PY = bar(distY; label="", title="distribution of Y", c=2)
PZ = bar(distZ; label="", title="distribution of Z", c=3)
plot(PX, PY, PZ; layout=(3, 1), ylim=(-0.01, 0.7))

# %%
distX = Categorical(0, 1/2, 1/6, 0, 1/3)
distY = Categorical(1/3, 0, 1/6, 1/2, 0)
distZ = Categorical(0, 1/6, 2/3, 1/6, 0)

@show probability_that_distX_wins_distY(distX, distY)
@show probability_that_distX_wins_distY(distY, distZ)
@show probability_that_distX_wins_distY(distZ, distX)
println()

PX = bar(distX; label="", title="distribution of X", c=1)
PY = bar(distY; label="", title="distribution of Y", c=2)
PZ = bar(distZ; label="", title="distribution of Z", c=3)
plot(PX, PY, PZ; layout=(3, 1), ylim=(-0.01, 0.7))

# %%
