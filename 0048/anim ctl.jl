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
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

function anim_lln_clt_rw(dist; nstart=1, nstep=5, nmax=1001, L=200,
        gn="lln_clt_rw.gif")
    μ, σ = mean(dist), std(dist)
    X = cumsum(rand(dist - μ, nmax, L); dims=1)
    anim = @animate for n in [nstart:nstep:nmax; fill(nmax, 30)]
        P = @views plot(1:n, X[1:n, :]; label="", alpha=0.5, lw=0.2)
        plot!(x-> 2σ*√x, 0, n; label="", c=:red, ls=:dot)
        plot!(x->-2σ*√x, 0, n; label="", c=:red, ls=:dot)
        Q = @views plot(1:n, X[1:n, :] ./ .√(1:n); label="", alpha=0.5, lw=0.2)
        plot!(x-> 2σ, 0, n; label="", c=:red, ls=:dot)
        plot!(x->-2σ, 0, n; label="", c=:red, ls=:dot)
        R = @views plot(1:n, X[1:n, :] ./ (1:n); label="", alpha=0.5, lw=0.2)
        plot!(x-> 2σ/√x, 1, n; label="", c=:red, ls=:dot)
        plot!(x->-2σ/√x, 1, n; label="", c=:red, ls=:dot)
        plot(P, Q, R; size=(500, 800), layout=(3, 1))
    end
    gif(anim, gn)
end

# %%
anim_lln_clt_rw(Exponential())

# %%
