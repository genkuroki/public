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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Distributions
using Optim
using Random
Random.seed!(4649373)
using StatsFuns
using StatsPlots
default(fmt=:png, titlefontsize=10, size=(400, 250))

logisticmodel(x, a, b) = Bernoulli(logistic(a + b*x))

randlogistic(x, a, b) = @. rand(logisticmodel(x, a, b))

neglogliklogistic(x, y, a, b) =
    -sum(logpdf(logisticmodel(x, a, b), y) for (x, y) in zip(x, y))

function logistic_regression(x, y; alg=LBFGS())
    f(w) = neglogliklogistic(x, y, w[1], w[2])
    o = optimize(f, zeros(2), alg)
    o.minimizer
end

function plot_logistic_regression(; a₀=1.0, b₀=2.0, x=rand(Uniform(-3, 3), 80), alg=LBFGS())
    y = randlogistic(x, a₀, b₀)

    @show a₀, b₀
    @show â, b̂ = logistic_regression(x, y; alg)

    a = range(-1, 3, 100)
    b = range(-1, 5, 100)
    negloglik = neglogliklogistic.(Ref(x), Ref(y), a', b)
    
    xlim = extrema(x)
    xlim = xlim[1] - 0.05*(xlim[2] - xlim[1]), xlim[2] + 0.05*(xlim[2] - xlim[1])
    P = scatter(x, y; ms=4, msc=:auto, alpha=0.7, label="data", legend=:bottomright)
    plot!(x -> logistic(â + b̂*x), xlim...; label="estimate", c=:red, lw=1.5)
    plot!(x -> logistic(a₀ + b₀*x), xlim...; label="true", ls=:dash, c=2)
    plot!(xguide="\$x\$", yguide="\$y\$")
    plot!(ytick=0:0.1:1)
    plot!(bottommargin=4Plots.mm)

    Q = heatmap(a, b, exp.(-negloglik); colorbar=false)
    scatter!([â], [b̂]; label="", c=1, msc=:auto)
    scatter!([a₀], [b₀]; label="", c=:lightblue, msc=:auto, marker=:star)
    plot!(xlim=extrema(a), ylim=extrema(b))
    plot!(xguide="\$a\$", yguide="\$b\$")

    plot(P, Q; size=(800, 250))
end

plot_logistic_regression()

# %%
plot_logistic_regression()

# %%
plot_logistic_regression()

# %%
plot_logistic_regression()

# %%
plot_logistic_regression()

# %%
plot_logistic_regression()

# %%
using Distributions
using Optim
using Random
Random.seed!(4649373)
using StatsFuns
using StatsPlots
default(fmt=:png, titlefontsize=10, size=(400, 250))

logisticmodel_2d(x1, x2, a, b1, b2) = Bernoulli(logistic(a + b1*x1 + b2*x2))

randlogistic_2d(x1, x2, a, b1, b2) = @. rand(logisticmodel_2d(x1, x2, a, b1, b2))

neglogliklogistic_2d(x1, x2, y, a, b1, b2) =
    -sum(logpdf(logisticmodel_2d(x1, x2, a, b1, b2), y) for (x1, x2, y) in zip(x1, x2, y))

function logistic_regression_2d(x1, x2, y; alg=LBFGS())
    f(w) = neglogliklogistic_2d(x1, x2, y, w[1], w[2], w[3])
    o = optimize(f, zeros(3), alg)
    o.minimizer
end

function logistic_regression_2d(; a0 = 1.0, b10 = 2.0, b20 = -2.0, n = 80, alg = LBFGS())
    x1, x2 = rand(Uniform(-2, 2), n), rand(Uniform(-2, 2), n)
    y = randlogistic_2d(x1, x2, a0, b10, b20)

    @show a0, b10, b20
    @show â, b̂₁, b̂₂ = logistic_regression_2d(x1, x2, y; alg)

    scatter(x1, x2; label="", c=:bwr, msc=:auto, marker_z=y, colorbar=false)
    plot!(x1 -> -(â + b̂₁*x1)/b̂₂; label="estimate", c=:red, lw=2)
    plot!(x1 -> -(a0 + b10*x1)/b20; label="true", c=2, ls=:dash)
    plot!(xlim=(-2, 2), ylim=(-2, 2), legend=:outertop)
    plot!(size=(300, 320))
end

logistic_regression_2d()

# %%
logistic_regression_2d()

# %%
logistic_regression_2d()

# %%
logistic_regression_2d()

# %%
logistic_regression_2d()

# %%
logistic_regression_2d()

# %%
