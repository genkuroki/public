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

# %% [markdown]
# See Section 2 of https://arxiv.org/abs/1203.3503
#
# <img src="2023-02-01.png">

# %% [markdown]
# <img src="IMG_8288.PNG" width=80%>

# %% [markdown]
# <img src="IMG_8289.PNG" width=80%>

# %%
using LinearAlgebra
using Distributions
using StatsPlots
default(fmt=:png)

# %%
function rand_XYZU(; c1, c0, c2)
    # Set var(U) = var(Z) = var(X) = 1.
    c3 = √(1 - c1^2)
    U = randn()
    Z = randn()
    X = c1*U + c3*Z
    Y = c0*X + c2*U
    [X, Y, Z, U]
end

function rand_XYZU(n; c1=0.5, c0=1, c2=1)
    samples = [rand_XYZU(; c1, c0, c2) for _ in 1:n]
    stack(samples)
end

function calc_all(; c1=0.5, c0=1, c2=1, λ=0, N=10^6)
    @show N
    c3 = √(1 - c1^2)
    @show c1 c3 c0 c2
    println()

    XYZU = rand_XYZU(N; c1, c0, c2)
    X, Y, Z, U = eachrow(XYZU)
    
    @show c0
    println()
    @show c0 + c1*c2
    @show a = [X;;] \ Y
    println()
    @show c0 + c1*c2/(1 - c3^2), -c1*c2*c3/(1 - c3^2)
    @show b = [X Z] \ Y
    println()
    @show c0*c1 + c2, c0*c3
    @show c = [U Z] \ Y
    
    nothing
end

# %% [markdown]
# <img src="2023-02-01_fig1.png" width="30%">

# %%
for c1 in -0.9:0.1:0.9
    calc_all(; c1)
    println("="^80)
end

# %% [markdown]
# <img src="2023-02-01_fig1.png" width="30%">

# %%
for c1 in -0.9:0.1:0.9
    calc_all(; c1, c2=0.2/c1)
    println("="^80)
end

# %% [markdown]
# <img src="2023-02-01_fig1.png" width="30%">

# %%
for t in range(-π, 0, 20)[begin+1:end-1]
    c1 = cos(t)
    calc_all(; c1, c2=0.2/c1)
    println("="^80)
end

# %% [markdown]
# <img src="2023-02-01_fig1.png" width="30%">

# %%
