# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://discourse.julialang.org/t/is-there-a-maximum-f-op-itrs-in-julia/63868/11

# %%
versioninfo()

# %%
using BenchmarkTools

function max_abs(a, b)
    m = zero(promote_type(eltype(a), eltype(b)))
    for i in eachindex(a, b)
        tmp = abs(a[i] - b[i]) 
        tmp > m && (m = tmp)
    end 
    m
end

N = 500
a = randn(N, N)
b = randn(N, N)

@show VERSION
@show Threads.nthreads()
println()
@show N
print("simple for loop:     ")
@btime max_abs($a, $b)
print("maximum(generator):  ")
@btime maximum(abs(i - j) for (i, j) in zip($a, $b))
print("maximum(abs, a - b): ")
@btime maximum(abs, $a - $b)
print("maximum splat(abs∘-):")
@btime maximum(Base.splat(abs∘-), zip($a, $b))
print("mapreduce abs∘- max: ")
@btime mapreduce(abs∘-, max, $a, $b)
using Tullio
print("Tullio:              ")
@btime @tullio (max) _ := abs($a[i] - $b[i])
print("Tullio (LoopVect.):  ")
using LoopVectorization
@btime @tullio (max) _ := abs($a[i] - $b[i])

function max_abs_turbo(a, b)
    m = zero(promote_type(eltype(a), eltype(b)))
    @turbo for i in eachindex(a, b)
        m = max(m, abs(a[i] - b[i]))
    end 
    m
end

function max_abs_tturbo(a, b)
    m = zero(promote_type(eltype(a), eltype(b)))
    @tturbo for i in eachindex(a, b)
        m = max(m, abs(a[i] - b[i]))
    end 
    m
end

print("LoopVect. @turbo:    ")
@btime max_abs_turbo($a, $b)
print("LoopVect. @tturbo:   ")
@btime max_abs_tturbo($a, $b);

# %%
max_abs(a, b) ==
maximum(abs(i - j) for (i, j) in zip(a, b)) ==
maximum(abs, a - b) ==
maximum(Base.splat(abs∘-), zip(a, b)) ==
mapreduce(abs∘-, max, a, b) ==
(@tullio (max) _ := abs(a[i] - b[i])) ==
max_abs_turbo(a, b) ==
max_abs_tturbo(a, b)

# %%
print("Tullio avx=false:    ")
@btime @tullio avx=false (max) _ := abs($a[i] - $b[i])
print("Tullio threads=false:")
@btime @tullio threads=false (max) _ := abs($a[i] - $b[i]);

# %%
