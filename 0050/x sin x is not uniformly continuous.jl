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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Plots
default(fmt=:png)

f(x) = x * sin(x)
suff_large_n(δ) = ceil(Int, (2/min(δ, π/2) + 1)/4)
plot(f, 0, 20π; label="x sin x", size=(400, 250)) |> display

PP = []
for δ in (4, 1, 0.5, 0.1)
    n = suff_large_n(δ)
    a = min(δ, π/2)
    P = plot(f, 2π*n - 1.05δ, 2π*n + 1.05δ; label=("x sin x"))
    plot!(f, 2π*n-a, 2π*n-a/2; label="", c=:blue, lw=3)
    plot!(f, 2π*n+a/2, 2π*n+a; label="", c=:blue, lw=3)
    hline!([0]; label="", c=2, ls=:dot)
    hline!([1, -1]; label="±1", c=2, ls=:dash)
    vline!([2π*n], label="2πn", c=3, ls=:dot)
    vline!([2π*n - a, 2π*n + a], label="2πn ± min(δ, π/2)", c=3, ls=:dash)
    vline!([2π*n - a/2, 2π*n + a/2], label="2πn ± min(δ, π/2)/2", c=3, ls=:dashdot)
    title!("δ = $δ,  n = $n")
    push!(PP, P)
end
plot(PP...; size=(800, 500), layout=(2, 2))

# %%
using Plots
default(fmt=:png)

a(n) = (-1)^(n-1) * (1 - 1/n)

N = 100
scatter(a.(n), zeros(length(n)); label="", msc=:auto, ms=3, alpha=0.7)
plot!(xlim=(-1.1, 1.1), ylim=(-1, 1), size=(600, 100), yaxis=false, ytick=false)
title!("N = $N") |> display

Ns = [vcat((fill(n, ceil(Int, 20/n)) for n in 1:100)...); fill(100, 20)]
anim = @animate for N in Ns
    n = 1:N
    scatter(a.(n), zeros(length(n)); label="", msc=:auto, ms=3, alpha=0.7)
    plot!(xlim=(-1.1, 1.1), ylim=(-1, 1), size=(600, 100), yaxis=false, ytick=false)
    title!("N = $N")
end
gif(anim, "aₙ.gif")

# %%
