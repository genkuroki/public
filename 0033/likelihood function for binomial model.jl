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
using StatsPlots
default(fmt=:png, titlefontsize=10)
plot(p -> p*(1-p), 0, 1; label="p -> p*(1-p)")

# %%
n, k = 100, 90
plot(legend=:topleft)
plot!(p -> pdf(Binomial(n, p), k), 0.7, 1;
    label="p -> pdf(Binomial(n, p), k)")
plot!(p -> pdf(Normal(n*p, √(n*p*(1-p))), k); 
    label="p -> pdf(Normal(n*p, √(n*p*(1-p))), k)", ls=:dash)

# %%
n, k = 1000, 990
plot(legend=:topleft)
plot!(p -> pdf(Binomial(n, p), k), 0.97, 1;
    label="p -> pdf(Binomial(n, p), k)")
plot!(p -> pdf(Normal(n*p, √(n*p*(1-p))), k);
    label="p -> pdf(Normal(n*p, √(n*p*(1-p))), k)", ls=:dash)

# %%
n, p = 100, 0.9
plot(legend=:topleft)
plot!(k -> pdf(Binomial(n, p), round(Int, k)), 70, 100;
    label="k -> pdf(Binomial(n, p), round(Int, k))")
plot!(k -> pdf(Normal(n*p, √(n*p*(1-p))), k);
    label="k -> pdf(Normal(n*p, √(n*p*(1-p))))", ls=:dash)

# %%
using Distributions
using StatsFuns
using StatsPlots
default(fmt=:png, titlefontsize=10)

n, k = 100, 90
plot(legend=:topleft)
plot!(x -> (p=logistic(x); pdf(Binomial(n, p), k)), logit(0.66), logit(0.985); label="")
plot!(x -> (p=logistic(x); p̂=k/n; pdf(Normal(x, 1/√(n*p*(1-p))), logit(p̂))*p̂*(1-p̂)); label="", ls=:dash)
plot!(x -> (p=logistic(x); p̂=k/n; pdf(Normal(x, 1/√(n*p*(1-p))), logit(p̂))*p*(1-p)); label="", ls=:dashdot)

# %%
using Distributions
using StatsFuns
using StatsPlots
default(fmt=:png, titlefontsize=10)

n, k = 100, 95
plot(legend=:topleft)
plot!(x -> (p=logistic(x); pdf(Binomial(n, p), k)), logit(0.7), logit(0.993); label="")

# %%
