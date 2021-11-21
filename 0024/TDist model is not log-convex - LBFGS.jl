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
#     display_name: Julia 1.6.4
#     language: julia
#     name: julia-1.6
# ---

# %%
using Distributions
using Plots
using Optim

traceon = Optim.Options(
    store_trace = true,
    extended_trace = true
)

# model
Distributions.TDist(μ, ρ, ν) = LocationScale(μ, ρ, TDist(ν))

# test data
X = [-0.01, 0.01, 1.0]

# estimations with different initial values
r1 = optimize(x -> -loglikelihood(TDist(x[1], eps() + 10^x[2], eps() + 10^x[3]), X), [0.33, 0.0, 0.16], LBFGS(), traceon)
r2 = optimize(x -> -loglikelihood(TDist(x[1], eps() + 10^x[2], eps() + 10^x[3]), X), [0.33, 0.0, 0.18], LBFGS(), traceon)
r3 = optimize(x -> -loglikelihood(TDist(x[1], eps() + 10^x[2], eps() + 10^x[3]), X), [0.33, 0.0, 0.27], LBFGS(), traceon);

# %%
function plot_result(r, X) 
    x = r.initial_x
    m = r.minimizer
    l = -minimum(r)
    t = "initial = $x\nfinal = $(round.(m; digits=3))\nfinal loglikelihood = $(round(l; digits=3))"
    xtick = axes(X, 1)
    normal = fit_mle(Normal, X)

    plot(; legend=:topleft)
    scatter!(xtick, X; label="data", xtick, xlim=(first(xtick)-0.5, last(xtick)+0.5))
    hline!([m[1]]; label="result of t-dist. model", lw=1.5)
    hline!([mean(normal)]; label="result of normal dist. model", ls=:dash, lw=1.5)
    title!(t; titlefontsize=10)
end

R1 = plot_result(r1, X)
R2 = plot_result(r2, X)
R3 = plot_result(r3, X)
plot(R1, R2, R3; size=(800, 700), xtick=1:3)

# %%
function plot_trace(r)
    x = r.initial_x
    c = hcat(x, (t.metadata["x"] for t in r.trace)...)
    l = -minimum(r)
    m = r.minimizer
    t = "initial = $x\nfinal = $(round.(m; digits=3))\nfinal loglikelihood = $(round(l; digits=3))"

    plot(; legend=:topleft)
    scatter3d!([x[1]], [x[2]], [x[3]]; label="initial", c=:blue, ms=3)
    plot!(c[1,:], c[2,:], c[3,:]; label="path", c=:darkcyan)
    scatter3d!([m[1]], [m[2]], [m[3]]; label="final", c=:yellow, ms=3)
    plot!(; xlabel="μ", ylabel="log₁₀ρ", zlabel="log₁₀ν")
    title!(t; titlefontsize=10)
end

P1 = plot_trace(r1)
P2 = plot_trace(r2)
P3 = plot_trace(r3)
plot(P1, P2, P3; size=(800, 700), tickfontsize=6)

# %%
function f(μ, log10ν, X)
    log10rhos = range(-15, 10; length=200)
    -maximum(loglikelihood(TDist(μ, 10^log10ρ, 10^log10ν), X) for log10ρ in log10rhos)
end

function plot_trace2d(r)
    x = r.initial_x
    c = hcat(x, (t.metadata["x"] for t in r.trace)...)
    l = -minimum(r)
    m = r.minimizer
    t = "initial = $x\nfinal = $(round.(m; digits=3))\nfinal loglikelihood = $(round(l; digits=3))"

    mus = range(-0.1, 0.4; length=200)
    log10nus = range(-2.5, 8.0; length=200)
    z = f.(mus', log10nus, Ref(X))

    plot(; legend=:topleft, colorbar=false)
    plot!(; xlim=extrema(mus), ylim=extrema(log10nus))
    heatmap!(mus, log10nus, z; clim=(-1, 5), c=reverse(cgrad(:CMRmap)))
    scatter!([x[1]], [x[3]]; label="initial", ms=4, c=:blue)
    plot!(c[1,:], c[3,:]; label="path", c=:cyan, lw=1.5)
    scatter!([m[1]], [m[3]]; label="final", ms=4, c=:yellow)
    plot!(; xlabel="μ", ylabel="log₁₀ν")
    title!(t; titlefontsize=10)
end

Q1 = plot_trace2d(r1)
Q2 = plot_trace2d(r2)
Q3 = plot_trace2d(r3)
plot(Q1, Q2, Q3; size=(800, 700))

# %%
