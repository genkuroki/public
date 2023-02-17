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
#     display_name: Julia 1.9.0-beta3
#     language: julia
#     name: julia-1.9
# ---

# %%
using LinearAlgebra
dot2(x) = dot(x, x)
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6)

# %%
n = 200
x1 = rand(Uniform(-√3, √3), n)
x2 = 10rand(Uniform(-√3, √3), n)
y = @. 1 + x1 + x2 + randn()

X = [ones(n) x1 x2]
@show β̂0, β̂1, β̂2 = X \ y
ŷ = @. β̂0 + β̂1*x1 + β̂2*x2
@show ŝ² = dot2(y - ŷ)/(n-3)

P1 = scatter(x1, y; label="(x1, y)", msc=:auto, alpha=0.5, ms=3)
P2 = scatter(x2, y; label="(x2, y)", msc=:auto, alpha=0.5, ms=3)

Q1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
Q2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(P1, P2, Q1, Q2; size=(800, 800), layout=(2, 2))

# %%
n = 200
x1 = rand(Uniform(-√3, √3), n)
x2 = 10rand(Uniform(-√3, √3), n)
y = @. 1 + x1 + x1^2 + x2 + randn()

X = [ones(n) x1 x2]
@show β̂0, β̂1, β̂2 = X \ y
ŷ = @. β̂0 + β̂1*x1 + β̂2*x2
@show ŝ² = dot2(y - ŷ)/(n-3)

P1 = scatter(x1, y; label="(x1, y)", msc=:auto, alpha=0.5, ms=3)
P2 = scatter(x2, y; label="(x2, y)", msc=:auto, alpha=0.5, ms=3)

Q1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
Q2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(P1, P2, Q1, Q2; size=(800, 800), layout=(2, 2))

# %%
X = [ones(n) x1 (x1 .^2) x2]
@show γ̂0, γ̂11, γ̂12, γ̂2 = X \ y
ŷ = @. γ̂0 + γ̂12*x1 + γ̂12*x1^2 + γ̂2*x2
@show ŝ² = dot2(y - ŷ)/(n-4)

Q1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
Q2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(P1, P2, Q1, Q2; size=(800, 800), layout=(2, 2))

# %%
