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
#     display_name: Julia 1.7.0-beta2
#     language: julia
#     name: julia-1.7
# ---

# %%
@show VERSION
using BenchmarkTools
pimc_sum(n) = 4sum(rand()^2 + rand()^2 ≤ 1 for _ in 1:n)/n
@btime pimc_sum(10^7)

# %%
@show VERSION
using BenchmarkTools, Random
pimc_sum(n, rng) = 4sum(rand(rng)^2 + rand(rng)^2 ≤ 1 for _ in 1:n)/n
@btime pimc_sum(10^7, $(Random.default_rng()))

# %%
@show VERSION
using BenchmarkTools, Random
function pimc_for(n, rng)
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end
@btime pimc_for(10^7, $(Random.default_rng()))

# %%
@show VERSION
@btime pimc_sum(10^7, $(MersenneTwister()))
@btime pimc_for(10^7, $(MersenneTwister()))

# %%
@show VERSION
f3(n) = sum(k^3 for k in 1:n)
f4(n) = sum(k^4 for k in 1:n)
@btime f3(10^6)
@btime f4(10^6)

# %%
@show VERSION
f3(n) = sum(k^3 for k in 1:n)
f4(n) = sum(k^4 for k in 1:n)
@btime f3($(Int128(10^6)))
@btime f4($(Int128(10^6)))

# %%
@code_typed f3(10^6)

# %%
@code_typed f4(10^6)

# %%

# %%
@show VERSION
using BenchmarkTools
pimc_sum(n) = 4sum(rand()^2 + rand()^2 ≤ 1 for _ in 1:n)/n
@btime pimc_sum(10^7)

# %%
@show VERSION
using BenchmarkTools, Random
pimc_sum(n, rng) = 4sum(rand(rng)^2 + rand(rng)^2 ≤ 1 for _ in 1:n)/n
@btime pimc_sum(10^7, $(Random.default_rng()))

# %%
@show VERSION
using BenchmarkTools, Random
function pimc_for(n, rng)
    c = 0
    for _ in 1:n
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/n
end
@btime pimc_for(10^7, $(Random.default_rng()))

# %%
@show VERSION
@btime pimc_sum(10^7, $(MersenneTwister()))
@btime pimc_for(10^7, $(MersenneTwister()))

# %%
@show VERSION
f3(n) = sum(k^3 for k in 1:n)
f4(n) = sum(k^4 for k in 1:n)
@btime f3(10^6)
@btime f4(10^6)

# %%
@show VERSION
f3(n) = sum(k^3 for k in 1:n)
f4(n) = sum(k^4 for k in 1:n)
@btime f3($(Int128(10^6)))
@btime f4($(Int128(10^6)))

# %%
@code_typed f3(10^6)

# %%
@code_typed f4(10^6)

# %%
