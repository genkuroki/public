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
default(fmt=:png, titlefontsize=8, tickfontsize=6, legendfontsize=7, guidefontsize=7)

# %%
function rand_XYZU(; c1, c0, c2, δ, ε)
    # Set var(U) = var(Z) = var(X) = 1.
    c3 = √(1 - (c1^2 + δ^2))
    U = randn()
    Z = randn()
    X = c1*U + c3*Z + δ*randn()
    Y = c0*X + c2*U + ε*randn()
    [X, Y, Z, U]
end

function rand_XYZU(n; c1=0.5, c0=1, c2=1, δ=0.2, ε=0.2)
    samples = [rand_XYZU(; c1, c0, c2, δ, ε) for _ in 1:n]
    stack(samples)
end

function calc_all(; c1=0.5, c0=1, c2=1, δ=0.1, ε=0.1, n=10^6)
    @show n
    c3 = √(1 - (c1^2 + δ^2))
    @show c0 c1 c2 c3 δ ε
    println()

    XYZU = rand_XYZU(n; c1, c0, c2, δ, ε)
    X, Y, Z, U = eachrow(XYZU)
    
    @show c0 + c1*c2
    @show a = [X;;] \ Y
    println()
    @show c0 + c1*c2/(1 - c3^2), -c1*c2*c3/(1 - c3^2)
    @show b = [X Z] \ Y
    println()
    @show c0*c1 + c2, c0*c3
    @show c = [U Z] \ Y
    println()
    @show c0, c2
    @show d = [X U] \ Y
    
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
# 誤差の比較

function calc_all(; c1=0.5, c0=1, c2=1, δ=0.2, ε=0.2, n=10^3, L=10^6, bin=200)
    @show L
    @show n
    c3 = √(1 - c1^2)
    @show c0 c1 c2 c3 δ ε
    println()

    XYZU = rand_XYZU(n; c1, c0, c2, δ, ε)
    X, Y, Z, U = eachrow(XYZU)
    
    @show c0 + c1*c2
    @show a = [X;;] \ Y
    println()
    @show c0 + c1*c2/(1 - c3^2), -c1*c2*c3/(1 - c3^2)
    @show b = [X Z] \ Y
    println()
    @show c0*c1 + c2, c0*c3
    @show c = [U Z] \ Y
    println()
    @show c0, c2
    @show d = [X U] \ Y
    
    XYZU_test = rand_XYZU(L; c1, c0, c2, δ, ε)
    error_a = [y - a[1]*x for (x, y, z, u) in eachcol(XYZU_test)]
    error_b = [y - (b[1]*x + b[2]*z) for (x, y, z, u) in eachcol(XYZU_test)]
    error_c = [y - (c[1]*u + c[2]*z) for (x, y, z, u) in eachcol(XYZU_test)]
    error_d = [y - (d[1]*x + d[2]*u) for (x, y, z, u) in eachcol(XYZU_test)]
    
    P1 = plot()
    stephist!(error_a; norm=true, bin, label="Y ~ X", c=1)
    stephist!(error_b; norm=true, bin, label="Y ~ X + Z", ls=:dash, c=2)
    title!("distribution of errors (sample size: n = $n)")
    
    P2 = plot()
    stephist!(error_a; norm=true, bin, label="Y ~ X", c=1)
    stephist!(error_b; norm=true, bin, label="Y ~ X + Z", ls=:dash, c=2)
    stephist!(error_c; norm=true, bin, label="Y ~ U + Z", ls=:dot, c=3)
    stephist!(error_d; norm=true, bin, label="Y ~ X + U", ls=:dashdot, c=4)
    title!("distribution of errors (sample size: n = $n)")
    
    plot(P1, P2; size=(900, 280))
end

# %% [markdown]
# <img src="2023-02-01_fig1.png" width="30%">

# %%
c1 = 0.05
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.1
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.2
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.3
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.4
calc_all(; c1, c2=0.2/c1, n=10^3)

# %% [markdown]
# <img src="2023-02-01_fig1.png" width="30%">

# %%
c1 = 0.5
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.6
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.7
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.8
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.9
calc_all(; c1, c2=0.2/c1, n=10^3)

# %%
c1 = 0.95
calc_all(; c1, c2=0.2/c1, n=10^3)

# %% [markdown]
# <img src="2023-02-01_fig1.png" width="30%">

# %%
