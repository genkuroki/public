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
# $$
# \begin{matrix}
# U          &          &     &          & V \\
# \downarrow & \searrow &     & \swarrow & \downarrow \\
# \downarrow &          & Z   &          & \downarrow \\
# \downarrow &          &     &          & \downarrow \\
# X          & \to      & \to & \to      & Y \\
# \end{matrix}
# \quad
# \begin{cases}
# Z = qU + sV + \varepsilon_1 \\
# X = pU + \varepsilon_2 \\
# Y = aX + rV + \varepsilon_3 \\
# \end{cases}
# \quad
# \begin{cases}
# \text{$U,V,\varepsilon_1,\varepsilon_2,\varepsilon_3$ are i.i.d $\sim$ Normal(0,1).} \\
# E[Z^2] = q^2 + s^2 + \sigma_1^2 := 1 \\
# E[X^2] = p^2 + \sigma_2^2 := 1 \\
# \end{cases}
# $$
#
# <img src="IMG_8351.PNG" width=80%>

# %%
using Distributions
using StatsPlots
default(fmt=:png, size=(400, 250),
    titlefontsize=8, tickfontsize=6, legendfontsize=7, guidefontsize=7)
rd(x) = round(x; digits=3)

# %%
function rand_XYZUV(; σ₁, σ₂, σ₃, a, p = √(1 - σ₂^2), q, r, s=√(1-q^2-σ₁^2))
    U = randn()
    V = randn()
    Z = q*U + s*V + σ₁*randn()
    X = p*U + σ₂*randn()
    Y = a*X + r*V + σ₃*randn()
    [X, Y, Z, U, V]
end

function rand_XYZUV(n;
        σ₁ = 0.2, σ₂ = 0.2, σ₃ = 0.2, 
        a = 1.0, p = √(1 - σ₂^2), q = √((1 - σ₁^2)/2), r = √10, s = √(1 - q^2 - σ₁^2))
    [rand_XYZUV(; σ₁, σ₂, σ₃, a, p, q, r) for _ in 1:n] |> stack
end

rand_XYZUV(100) .|> rd

# %%
function show_XYZUV(n = 10^6;
        σ₁ = 0.2, σ₂ = 0.2, σ₃ = 0.2, 
        a = 1.0, p = √(1 - σ₂^2), q = √((1 - σ₁^2)/2), r = √10, s = √(1 - q^2 - σ₁^2))
    println("====== model parameters")
    @show rd(σ₁) rd(σ₂) rd(σ₃) rd(a) rd(p) rd(q) rd(r) rd(s) n
    println()

    data = rand_XYZUV(n; σ₁, σ₂, σ₃, a, p, q, r, s)
    X, Y, Z, U, V = eachrow(data)
    println("====== summary of data")
    @show rd.(mean.((X, Y, Z, U, V)))
    @show rd.(var.((X, Y, Z, U, V)))
    @show rd.(std.((X, Y, Z, U, V)))
    println()

    println("====== resgressions")
    println("=== Y ~ α₁X")
    @show rd(a)
    α = [X;;] \ Y
    @show rd.(α)
    println("=== Y ∼ β₁X + β₂Z")
    @show rd(a - p*q*r*s/(1 - p^2*q^2)), rd(r*s/(1 - p^2*q^2))
    β = [X Z] \ Y
    @show rd.(β)
    println("=== Y ∼ γ₁X + γ₄V")
    @show rd(a), rd(r)
    γ = [X V] \ Y
    @show rd.(γ)
    println()
    
    println("====== M-bias")
    @show rd(β[1] - a)
    println()

    println("====== errors of prediction for Y")
    err_α = @. Y - α[1]*X
    err_β = @. Y - (β[1]*X + β[2]*Z)
    err_γ = @. Y - (γ[1]*X + γ[2]*V)
    @show rd.(mean.((err_α, err_β, err_γ)))
    @show rd.(var.((err_α, err_β, err_γ)))
    @show rd.(std.((err_α, err_β, err_γ)))
    println()

    P1 = plot()
    stephist!(err_α; norm=true, label="Y ∼ X", ls=:solid)
    stephist!(err_β; norm=true, label="Y ∼ X + Z", ls=:dash)
    plot!(xguide="error of prediction for Y")
    
    P2 = plot()
    stephist!(err_α; norm=true, label="Y ∼ X", ls=:solid)
    stephist!(err_β; norm=true, label="Y ∼ X + Z", ls=:dash)
    stephist!(err_γ; norm=true, label="Y ∼ X + V", ls=:dashdot)
    plot!(xguide="error of prediction for Y")

    plot(P1, P2; size=(800, 250))
    plot!(bottommargin=4Plots.mm)
end

# %%
show_XYZUV()

# %% [markdown]
# $$
# \begin{matrix}
# U          &          &     &          & V \\
# \downarrow & \searrow &     & \swarrow & \downarrow \\
# \downarrow &          & Z   &          & \downarrow \\
# \downarrow &          &     &          & \downarrow \\
# X          & \to      & \to & \to      & Y \\
# \end{matrix}
# \quad
# \begin{cases}
# Z = qU + sV + \varepsilon_1 \\
# X = pU + \varepsilon_2 \\
# Y = aX + rV + \varepsilon_3 \\
# \end{cases}
# \quad
# \begin{cases}
# \text{$U,V,\varepsilon_1,\varepsilon_2,\varepsilon_3$ are i.i.d $\sim$ Normal(0,1).} \\
# E[Z^2] = q^2 + s^2 + \sigma_1^2 := 1 \\
# E[X^2] = p^2 + \sigma_2^2 := 1 \\
# \end{cases}
# $$

# %%
