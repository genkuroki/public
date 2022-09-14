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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions
using Optim
using QuadGK
using Random
using StatsPlots
default(fmt=:png, titlefontsize=10, size=(500, 300))

# %% [markdown]
# $X,Y\sim\mathrm{Uniform}(0,1)$ が独立ならば $Z=X/Y$ が従う分布の密度函数は次になる:
#
# $$
# f(z) = \begin{cases}
# 0 & (z < 0) \\
# 1/2 & (0\le z\le 1) \\
# 1/(2z^2) & (1 \le z) \\
# \end{cases}
# $$

# %%
f(z) = z < 0 ? 0.0 : 0 ≤ z ≤ 1 ? 1/2 : 1/(2z^2)
quadgk(f, 0, 1, Inf)

# %%
n = 10^6
X = rand(Uniform(0, 1), n)
Y = rand(Uniform(0, 1), n)
Z = X ./ Y
stephist(Z; norm=true, bin=[0:0.1:30; Inf], label="X/Y")
plot!(f, -0.1, 30; label="f(z)", xlim=(-0.5, 10), ls=:dash)

# %% [markdown]
# $X,Y\sim\mathrm{Uniform}(-1,1)$ が独立ならば $Z=X/Y$ が従う分布の密度函数は次になる:
#
# $$
# g(z) = \begin{cases}
# 1/4 & (|z|\le 1) \\
# 1/(4z^2) & (1 \le |z|) \\
# \end{cases}
# $$

# %%
g(z) = abs(z) ≤ 1 ? 1/4 : 1/(4z^2)
quadgk(f, -Inf, -1, 0, 1, Inf)

# %%
n = 10^6
X = rand(Uniform(-1, 1), n)
Y = rand(Uniform(-1, 1), n)
Z = X ./ Y
stephist(Z; norm=true, bin=[-Inf; -30:0.1:30; Inf], label="Z = X/Y")
plot!(g, -30, 30; label="g(z)", xlim=(-10, 10), ls=:dash)
plot!(Cauchy(0, 0.8), -10, 10; label="Cauchy(0, 0.8)", ls=:dashdot)
plot!(xlim=(-10, 10))
title!("X, Y ~ Uniform(-1, 1), independently")

# %%
n, L = 10, 10^5
tmpX, tmpY = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
@time meanZ = [mean(((x, y),) -> x/y, zip(rand!(Uniform(-1, 1), tmpX), rand!(Uniform(-1, 1), tmpY))) for _ in 1:L]

stephist(meanZ; norm=true, bin=[-Inf; -10^4:0.1:10^4; Inf], label="mean(X/Y), n=$n")
plot!(Cauchy(0, 0.8), -10, 10; label="Cauchy(0, 0.8)", ls=:dashdot)
plot!(xlim=(-10, 10))
title!("X, Y ~ Uniform(-1, 1), independently")

# %%
n, L = 100, 10^5
tmpX, tmpY = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
@time meanZ = [mean(((x, y),) -> x/y, zip(rand!(Uniform(-1, 1), tmpX), rand!(Uniform(-1, 1), tmpY))) for _ in 1:L]

stephist(meanZ; norm=true, bin=[-Inf; -10^4:0.1:10^4; Inf], label="mean(X/Y), n=$n")
plot!(Cauchy(0, 0.8), -10, 10; label="Cauchy(0, 0.8)", ls=:dashdot)
plot!(xlim=(-10, 10))
title!("X, Y ~ Uniform(-1, 1), independently")

# %%
n, L = 1000, 10^5
tmpX, tmpY = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
@time meanZ = [mean(((x, y),) -> x/y, zip(rand!(Uniform(-1, 1), tmpX), rand!(Uniform(-1, 1), tmpY))) for _ in 1:L]

stephist(meanZ; norm=true, bin=[-Inf; -10^4:0.1:10^4; Inf], label="mean(X/Y), n=$n")
plot!(Cauchy(0, 0.8), -10, 10; label="Cauchy(0, 0.8)", ls=:dashdot)
plot!(xlim=(-10, 10))
title!("X, Y ~ Uniform(-1, 1), independently")

# %%
n, L = 10000, 10^5
tmpX, tmpY = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
@time meanZ = [mean(((x, y),) -> x/y, zip(rand!(Uniform(-1, 1), tmpX), rand!(Uniform(-1, 1), tmpY))) for _ in 1:L]

stephist(meanZ; norm=true, bin=[-Inf; -10^4:0.1:10^4; Inf], label="mean(X/Y), n=$n")
plot!(Cauchy(0, 0.8), -10, 10; label="Cauchy(0, 0.8)", ls=:dashdot)
plot!(xlim=(-10, 10))
title!("X, Y ~ Uniform(-1, 1), independently")

# %%
n = 10^6
X = rand(Uniform(-1, 1), n)
Y = rand(Uniform(-1, 1), n)
Z = X ./ Y

o = optimize(w -> -loglikelihood(Cauchy(w...), Z), [0.0, 1.0])
@show o
@show o.minimizer;

# %%
n, L = 10^4, 10^6
tmpX, tmpY = Vector{Float64}(undef, n), Vector{Float64}(undef, n)
@time meanZ = [mean(((x, y),) -> x/y, zip(rand!(Uniform(-1, 1), tmpX), rand!(Uniform(-1, 1), tmpY))) for _ in 1:L]

@time o = optimize(w -> -loglikelihood(Cauchy(w...), meanZ), [0.0, 1.0])
@show o
@show o.minimizer;

# %%
