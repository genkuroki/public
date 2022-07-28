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
default(fmt=:png, titlefontsize=12, legendfontsize=12)

# %%
n = 10^6
a, b, θ = 6, 14, 3
X = rand(Gamma(a, θ), n)
Y = rand(Gamma(b, θ), n)
P = @. X/(X+Y)
U = @. X/Y;

stephist(P; norm=true, label="X/(X+Y)")
plot!(Beta(a, b); label="Beta(a,b)", ls=:dash)
title!("X∼Gamma(a=$a, θ=$θ), Y∼Gamma(b=$b, θ=$θ)") |> display

stephist(U; norm=true, label="X/Y", xlim=(0, 2))
plot!(BetaPrime(a, b); label="BetaPrime(a,b)", ls=:dash)
title!("X∼Gamma(a=$a, θ=$θ), Y∼Gamma(b=$b, θ=$θ)") |> display

# %%
n = 10^6
ν₁, ν₂ = 12, 28
X = rand(Chisq(ν₁), n)
Y = rand(Chisq(ν₂), n)
F = @. (X/ν₁)/(Y/ν₂)

stephist(F; norm=true, label="X/Y", xlim=(0, 2*ν₂/ν₁))
plot!(FDist(ν₁, ν₂); label="FDist(ν₁, ν₂)", ls=:dash)
title!("X∼Chisq(ν₁=$ν₁), Y∼Chisq(ν₂=$ν₂)")

# %%
n = 10^6
ν = 3
Z = rand(Normal(0,1), n)
W = rand(Chisq(ν), n)
T = @. Z/√(W/ν)

stephist(T; norm=true, label="Z/√(W/ν)", xlim=(-6, 6))
plot!(TDist(ν); label="TDist(ν)", ls=:dash)
plot!(Normal(0,1), -6, 6; label="Normal(0,1)", ls=:dashdot)
title!("Z∼Normal(0,1), W∼Chisq(ν=$ν)")

# %%
