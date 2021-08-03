# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using Plots
x = rand(10^6)
b = 10.0 .^ (-5:0)
histogram(x; alpha=0.3, label="", bin=b, norm=true)

# %%
using Plots, StatsBase
x = rand(10^6)
b = 10.0 .^ (-5:0)
h = fit(Histogram{Float64}, x, b)
h.weights ./= diff(h.edges[1])
plot(h; alpha=0.3, label="")

# %%
using Plots
n = 10^6
f(x) = n * x * log(10)
x = rand(n)
b = 10.0 .^ (-5:0)
histogram(x; alpha=0.3, label="", bin=b, xscale=:log10, xlim=extrema(b))
x = 10.0 .^ (-5:0.01:0)
plot!(x, f; label="n × pdf on log-scaled xaxis", legend=:topleft)

# %%
using Plots
n = 10^6
f(x) = n * x * log(10)
x = rand(n)
b = 10.0 .^ (-5:0.2:0)
histogram(x; alpha=0.3, label="", bin=b, xscale=:log10, xlim=extrema(b))
x = 10.0 .^ (-5:0.01:0)
plot!(x, f; label="n × pdf on log-scaled xaxis", legend=:topleft)

# %%
using Plots
n = 10^6
f(x) = n * x * log(10)
x = rand(n)
b = 10.0 .^ (-5:0.2:0)
histogram(x; alpha=0.3, label="", bin=b, xscale=:log10, xlim=extrema(b))
x = 10.0 .^ (-5:0.01:0)
plot!(x, f; label="n × pdf on log-scaled xaxis", legend=:topleft)

# %%
using Plots, StatsBase
n = 10^6
f(x) = n * x * log(10)
x = rand(n)
b = 10.0 .^ (-5:0.2:0)
h = fit(Histogram{Float64}, x, b)
h.weights ./= diff(log10.(h.edges[1]))
plot(h; alpha=0.3, label="", xscale=:log10, xlim=extrema(b))
x = 10.0 .^ (-5:0.01:0)
plot!(x, f; label="n × pdf on log-scaled xaxis", legend=:topleft)

# %%
using Plots, StatsBase, Random
X = 3randexp(10^6)
bin = [0; [2^(k/2) - 2^(-1/2) for k in 0:8]]
@show round.(bin; digits=2)
h = fit(Histogram{Float64}, X, bin)
h.weights ./= diff(h.edges[1])
plot(h; alpha=0.3, label="")

# %%
using Plots, StatsBase, Random
n = 10^6
f(x) = n * exp(-x/3)/3 * x * log(10)
X = 3randexp(10^6)
bin = [1e-3; [2^(k/2) - 2^(-1/2) for k in 0:8]]
@show round.(bin; digits=2)
h = fit(Histogram{Float64}, X, bin)
h.weights ./= diff(log10.(h.edges[1]))
plot(h; alpha=0.3, label="", xscale=:log10, xlim=extrema(bin))
x = exp.(range(log.(extrema(bin))...; length=500))
plot!(x, f; label="n × pdf on log-scaled xaxis", legend=:topleft)

# %%
