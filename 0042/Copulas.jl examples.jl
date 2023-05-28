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

function plot_XY(X, Y; alpha=0.3, kwargs...)
    @show cor(X, Y)
    n = length(X)
    t = 1:min(n, 10^4)
    P1 = scatter(X[t], Y[t]; label="", ms=1.5, msc=:auto, alpha)
    plot!(xguide="X", yguide="Y")
    P2 = stephist(X; label="X", norm=true)
    plot!(Normal(); label="standard normal", ls=:dash)
    P3 = stephist(Y; label="Y", norm=true)
    plot!(Normal(); label="standard normal", ls=:dash)
    P4 = stephist(X+Y; label="X+Y", norm=true)
    plot!(fit(Normal, X+Y); label="normal", ls=:dash)
    plot!(; kwargs...)
    plot(P1, P3, P2, P4; layout=(2, 2), size=(800, 800))
end

function plot_XY(copula; n=10^6)
    mvdist = SklarDist(copula, (Normal(), Normal()))
    X, Y = eachrow(rand(mvdist, n))
    plot_XY(X, Y)
end

# %%
plot_XY(JoeCopula(2, 2.2))

# %%
plot_XY(TCopula(1, [1 0; 0 1]))

# %%
plot_XY(ClaytonCopula(2, 3.0))

# %%
RademacherDist = 2Bernoulli(1/2) - 1

n = 10^6
X = rand(Normal(), n)
Z = rand(RademacherDist, n)
Y = @. Z*X
plot_XY(X, Y; alpha=0.15, ylim = (-0.03, 1))

# %%
