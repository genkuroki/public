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
r(x) = round(x; digits=3)

# %% [markdown]
# https://arxiv.org/abs/1203.3503 Section 4
#
# $$
# \begin{matrix}
# Z          & & \\
# \downarrow & \searrow & \\
# X          & \to      & Y \\
# \uparrow   & \nearrow & \\
# U          & & \\
# \end{matrix}
# \qquad
# \begin{cases}
# X = bZ + pU + \varepsilon_1 \\
# Y = aX + cZ + qU + \varepsilon_2 \\
# \end{cases}
# \qquad
# E[X^2] = b^2+p^2+\sigma_1^2 = 1
# $$
#
# <img src="Section 4 of arxiv 1203.3503.jpg" width=80%>

# %%
function rand_XYZU(; σ₁=0.2, σ₂=0.2, a=1.0, b=0.9, c=1.0, p=√(1-b^2-σ₁^2), q=1.0)
    Z = randn()
    U = randn()
    X = b*Z + p*U + σ₁*randn()
    Y = a*X + c*Z + q*U + σ₂*randn()
    [X, Y, Z, U]
end

function rand_XYZU(n; σ₁=0.2, σ₂=0.2, a=1.0, b=0.9, c=0.5, p=√(1-b^2-σ₁^2), q=0.5)
    [rand_XYZU(; σ₁, σ₂, a, b, c, p, q) for _ in 1:n] |> stack
end

rand_XYZU(100) .|> r

# %%
function show_XYZU(; σ₁=0.2, σ₂=0.2, a=1.0, b=0.9, c=0.5, p=√(1-b^2-σ₁^2), q=0.5, n=10^6)
    @show r(σ₁) r(σ₂) r(a) r(b) r(c) r(p) r(q) n
    println()

    data = rand_XYZU(n; σ₁, σ₂, a, b, c, p, q)
    X, Y, Z, U = eachrow(data)
    @show mean.((X, Y, Z, U)) .|> r
    @show var.((X, Y, Z, U)) .|> r
    @show std.((X, Y, Z, U)) .|> r
    println()

    # Y ∼ α₁X
    @show a + b*c + p*q |> r
    α = [X;;] \ Y
    @show r.(α)
    println()
    # Y ∼ β₁X + β₂Z
    @show r(a + p*q/(1 - b^2)) r(c - b*p*q/(1 - b^2))
    β = [X Z] \ Y
    @show r.(β)
    println()
    # Y ∼ γ₁X + γ₂Z + γ₃U
    @show a, c, q
    γ = [X Z U] \ Y
    @show r.(γ)
    println()

    err_α = @. Y - α[1]*X
    err_β = @. Y - (β[1]*X + β[2]*Z)
    err_γ = @. Y - (γ[1]*X + γ[2]*Z + γ[3]*U)
    @show mean.((err_α, err_β, err_γ)) .|> r
    @show var.((err_α, err_β, err_γ)) .|> r
    @show std.((err_α, err_β, err_γ)) .|> r
    println()

    plot()
    stephist!(err_α; norm=true, label="Y ∼ X", ls=:solid)
    stephist!(err_β; norm=true, label="Y ∼ X + Z", ls=:dash)
    stephist!(err_γ; norm=true, label="Y ∼ X + Z + U", ls=:dashdot)
    plot!(xguide="error")
end

# %%
σ₁=0.2
σ₂=0.2
a=1.0
b=0.9
c=0.5
p=√(1-b^2-σ₁^2)
q=0.5
n = 10^6
show_XYZU(; σ₁, σ₂, a, b, c, p, q, n)

# %%
σ₁=0.2
σ₂=0.2
a=1.0
b=√((1 - σ₁^2)/2)
c=0.5
p=√(1-b^2-σ₁^2)
q=0.5
n = 10^6
show_XYZU(; σ₁, σ₂, a, b, c, p, q, n)

# %%
σ₁=0.2
σ₂=0.2
a=1.0
b=0.1
c=0.5
p=√(1-b^2-σ₁^2)
q=0.5
n = 10^6
show_XYZU(; σ₁, σ₂, a, b, c, p, q, n)

# %%
σ₁=0.4
σ₂=2.0
a=1.0
b=0.1
c=1/b
p=√(1-b^2-σ₁^2)
q=1/p
n = 10^6
show_XYZU(; σ₁, σ₂, a, b, c, p, q, n)

# %%
