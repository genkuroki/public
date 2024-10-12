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
#     display_name: Julia 1.11.0
#     language: julia
#     name: julia-1.11
# ---

# %%
using BenchmarkTools

soft_threashold_bc(x, λ) = @. sign(x) * max(0, abs(x) - λ)

function soft_threashold_naive(x, λ)
    y = zero(x)
    ind = findall(>(λ), x)
    y[ind] = x[ind] .- λ
    ind = findall(<(-λ),  x)
    y[ind] = x[ind] .+ λ
    y
end

x20 = randn(20)
@show soft_threashold_bc(x20, 1.0) == soft_threashold_naive(x20, 1.0)
print("soft_threashold_bc(x20, 1.0):    "); @btime soft_threashold_bc(x20, 1.0)
print("soft_threashold_naive(x20, 1.0) "); @btime soft_threashold_naive(x20, 1.0)
println()

x1000 = randn(1000)
@show soft_threashold_bc(x1000, 1.0) == soft_threashold_naive(x1000, 1.0)
print("soft_threashold_bc(x1000, 1.0):     "); @btime soft_threashold_bc(x1000, 1.0)
print("soft_threashold_naive(x1000, 1.0) "); @btime soft_threashold_naive(x1000, 1.0)
;

# %%
soft_threashold_bc!(y, x, λ) = @. y = sign(x) * max(0, abs(x) - λ)

y1000 = similar(x1000)
@show soft_threashold_bc(x1000, 1.0) == soft_threashold_bc!(y1000, x1000, 1.0)
print("soft_threashold_bc(x1000, 1.0):         "); @btime soft_threashold_bc($x1000, 1.0)
print("soft_threashold_bc!(y1000, x1000, 1.0): "); @btime soft_threashold_bc!($y1000, $x1000, 1.0);

# %%
soft_threashold_bc2!(y, x, λ) = @. y = (x > λ) * (x - λ) + (x < -λ) * (x + λ)

y1000 = similar(x1000)
z1000 = similar(x1000)
@show soft_threashold_bc!(y1000, x1000, 1.0) == soft_threashold_bc2!(z1000, x1000, 1.0)
print("soft_threashold_bc(x1000, 1.0):          "); @btime soft_threashold_bc($x1000, 1.0)
print("soft_threashold_bc!(y1000, x1000, 1.0):  "); @btime soft_threashold_bc!($y1000, $x1000, 1.0)
print("soft_threashold_bc2!(z1000, x1000, 1.0): "); @btime soft_threashold_bc2!($z1000, $x1000, 1.0);

# %%
@show findall(>(1.0), x1000) |> typeof
@show (x1000 .> 1.0) |> typeof
println()

function soft_threashold_naive2(x, λ)
    y = zero(x)
    ind = x .> λ
    y[ind] = x[ind] .- λ
    ind = x .< -λ
    y[ind] = x[ind] .+ λ
    y
end

@show soft_threashold_naive(x1000, 1.0) == soft_threashold_naive2(x1000, 1.0)
print("soft_threashold_naive(x1000, 1.0):  "); @btime soft_threashold_naive($x1000, 1.0)
print("soft_threashold_naive2(x1000, 1.0): "); @btime soft_threashold_naive2($x1000, 1.0)
;

# %%
versioninfo()

# %%
