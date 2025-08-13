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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
using Printf
using Plots
default(fmt=:png, size=(600, 260))

# %%
f(t, x, T, λ, θ₀, A) = A * sin(2π*(t/T - x/λ) + θ₀)

# %%
T = 1
λ = 1
θ₀ = 0
A = 1

@gif for t in range(0, 6, 201)
    plot(x -> f(t, x, T, λ, θ₀, A), -1, 3; label="")
    plot!(xguide="x", yguide="y")
    title!("T=$T, λ=$λ, t=" * @sprintf("%4.3f", t))
end

# %%
T = 1
λ = 2
θ₀ = 0
A = 1

@gif for t in range(0, 6, 201)
    plot(x -> f(t, x, T, λ, θ₀, A), -1, 3; label="")
    plot!(xguide="x", yguide="y")
    title!("T=$T, λ=$λ, t=" * @sprintf("%4.3f", t))
end

# %%
T = 1
λ = 1
θ₀ = 0
A = 1

@gif for t in range(0, 6, 201)
    plot(x -> f(t, x, T, λ, θ₀, A), -1, 3; label="")
    plot!(xguide="x", yguide="y")
    title!("T=$T, λ=$λ, t=" * @sprintf("%4.3f", t))
end

# %%
T = 2
λ = 1
θ₀ = 0
A = 1

@gif for t in range(0, 6, 201)
    plot(x -> f(t, x, T, λ, θ₀, A), -1, 3; label="")
    plot!(xguide="x", yguide="y")
    title!("T=$T, λ=$λ, t=" * @sprintf("%4.3f", t))
end

# %%
T = 1
λ = 1
θ₀ = 0
A = 1

@gif for t in range(0, 6, 201)
    plot(x -> f(t, x, T, λ, θ₀, A), -1, 3; label="")
    plot!(xguide="x", yguide="y")
    title!("T=$T, λ=$λ, t=" * @sprintf("%4.3f", t))
end

# %%
T = 2
λ = 2
θ₀ = 0
A = 1

@gif for t in range(0, 6, 201)
    plot(x -> f(t, x, T, λ, θ₀, A), -1, 3; label="")
    plot!(xguide="x", yguide="y")
    title!("T=$T, λ=$λ, t=" * @sprintf("%4.3f", t))
end

# %%
T = 1
λ = 1
θ₀ = 0
A = 1
t = range(0, 6, 201)
x = range(-1, 3, 201)
surface(t, x, (t, x) -> f(t, x, T, λ, θ₀, A); colorbar=false, camera=(30, 80))
plot!(xguide="t", yguide="x", zguide="y")
plot!(size=(600, 600))

# %%
T = 1
λ = 1
θ₀ = 0
A = 1
t = range(0, 6, 201)
x = range(-1, 3, 201)
heatmap(t, x, (t, x) -> f(t, x, T, λ, θ₀, A); colorbar=false, camera=(30, 80))
plot!(xguide="t", yguide="x")
plot!(size=(500, 500))

# %%
