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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions
using StatsPlots
using QuadGK
default(fmt=:png)

# %%
function plot_logFDist(m, n;
        fdist = FDist(m, n),
        xlim = log.(quantile.(fdist, (0.00001, 0.99999))))
    f(t) = pdf(fdist, exp(t))*exp(t)
    logf(t) = logpdf(fdist, exp(t)) + t
    μ = quadgk(t -> t*f(t), xlim...)[1]
    σ² = quadgk(t -> (t - μ)^2*f(t), xlim...)[1]
    @show normal = Normal(μ, √σ²)
    plot(f ; label="log FDist($m, $n)", xlim)
    plot!(normal; label="normal approx", ls=:dash)
    if n > 4
        @show μf, σf = mean(fdist), std(fdist)
        @show normal_delta = Normal(log(μf) - σf^2/(2μf), σf/μf)
        plot!(normal_delta; label="delta method", ls=:dashdot)
    end
    plot!()
end

# %%
plot_logFDist(3, 1.5)

# %%
plot_logFDist(10, 5)

# %%
plot_logFDist(30, 15)

# %%
plot_logFDist(100, 50)

# %%
plot_logFDist(300, 150)

# %%
