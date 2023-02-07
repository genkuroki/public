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
using Distributions
using StatsPlots
default(fmt=:png, size=(400, 250),
    titlefontsize=8, tickfontsize=6, legendfontsize=7, guidefontsize=7)
rd(x) = round(x; digits=3)

# %% [markdown]
# $$
# \begin{matrix}
#   &       &          & U   & & \\
#   &       & \swarrow &     & \searrow & \\
# Z & \to X & \to      & \to & \to & Y \\
# \end{matrix}
# \quad
# \begin{cases}
# X = bZ + pU + \sigma_1\varepsilon_1 \\
# Y = aX + qU + \sigma_2\varepsilon_2 \\
# \end{cases}
# \quad
# \begin{cases}
# \text{$\varepsilon_1,\varepsilon_2,Z$ are independent and have mean zero.} \\
# E[X^2] = b^2 + p^2 + \sigma_1^2 := 1 \\
# \end{cases}
# $$

# %%
function rand_XYZU(; σ₁, σ₂, a, b, p = √(1-b^2-σ₁^2), q)
    U = randn()
    Z = randn()
    X = b*Z + p*U + σ₁*randn()
    Y = a*X + q*U + σ₂*randn()
    [X, Y, Z, U]
end

function rand_XYZU(n;
        σ₁ = 0.2, σ₂ = 0.2, 
        a = 1.0, b = √((1 - σ₁^2)/2), p = √(1-b^2-σ₁^2), q = 1/p)
    [rand_XYZU(; σ₁, σ₂, a, b, p, q) for _ in 1:n] |> stack
end

rand_XYZU(100) .|> rd

# %%
function show_XYZU(n = 10^6;
        σ₁ = 0.2, σ₂ = 0.2, 
        a = 1.0, b = √((1 - σ₁^2)/2), p = √(1-b^2-σ₁^2), q = 1/p)
    println("====== model parameters")
    @show rd(σ₁) rd(σ₂) rd(a) rd(b) rd(p) rd(q) n
    println()

    data = rand_XYZU(n; σ₁, σ₂, a, b, p, q)
    X, Y, Z, U = eachrow(data)
    println("====== summary of data")
    @show rd.(mean.((X, Y, Z, U)))
    @show rd.(var.((X, Y, Z, U)))
    @show rd.(std.((X, Y, Z, U)))
    println()

    println("====== resgressions")
    println("=== Y ~ α₁X")
    α = [X;;] \ Y
    @show rd(a + p*q)
    @show rd.(α)
    println("=== Y ∼ β₁Z")
    β = [Z;;] \ Y
    @show rd.(a*b)
    @show rd.(β)
    println("=== X ∼ γ₁Z")
    γ = [Z;;] \ X
    @show rd.(b)
    @show rd.(γ)
    println("=== instrumental variable method")
    @show rd(a)
    @show rd(β[1]/γ[1])
    println()

    println("====== errors of prediction for Y")
    err_α = @. Y - α[1]*X
    err_β = @. Y - β[1]/γ[1]*Z
    @show rd.(mean.((err_α, err_β)))
    @show rd.(var.((err_α, err_β)))
    @show rd.(std.((err_α, err_β)))
    println()

    plot()
    stephist!(err_α; norm=true, label="Y ∼ X", ls=:solid)
    stephist!(err_β; norm=true, label="Y ∼ Z", ls=:dash)
    plot!(xguide="error of prediction for Y")
end

# %%
show_XYZU(a=2)

# %%
show_XYZU(a=2, b=0.4)

# %%
show_XYZU(a=2.0, b=0.97)

# %%
