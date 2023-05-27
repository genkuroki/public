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
using Copulas
using Random
using StatsPlots
default(fmt=:png)

function plot_XY(copula; n=10^6)
    mvdist = SklarDist(copula, (Normal(), Normal()))
    X, Y = eachrow(rand(mvdist, n))
    t = 1:min(n, 10^4)
    P1 = scatter(X[t], Y[t]; label="", ms=1.5, msc=:auto, alpha=0.5)
    P2 = stephist(X; label="X", norm=true)
    plot!(Normal(); label="", ls=:dash)
    P3 = stephist(Y; label="Y", norm=true)
    plot!(Normal(); label="", ls=:dash)
    P4 = stephist(X+Y; label="X+Y", norm=true)
    plot!(fit(Normal, X+Y); label="", ls=:dash)
    plot(P1, P3, P2, P4; layout=(2, 2), size=(800, 800))
end

# %%
plot_XY(JoeCopula(2, 2.2))

# %%
plot_XY(TCopula(2, [1 0; 0 1]))

# %%
plot_XY(ClaytonCopula(2, 3.0))

# %%
