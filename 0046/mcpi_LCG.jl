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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
module O

export rand_LCG, seed_LCG!

# Linear congruential generators (LCGs)
# Parameters are provided by Park and Miller
# See https://c-faq.com/lib/rand.html

mutable struct LinearCongruentialGenerator <: Function
    seed::Int32
end

const a = Int32(48271)
const m = Int32(2147483647)
const q = m ÷ a
const r = m % a

@inline function (f::LinearCongruentialGenerator)()
    seed = f.seed
    hi, lo = divrem(seed, q)
    seed = a * lo - r * hi
    seed = ifelse(seed > zero(seed), seed, seed + m)
    f.seed = seed
    seed / m
end

const rand_LCG = LinearCongruentialGenerator(Int32(20231226))

seed_LCG!(seed) = rand_LCG.seed = Int32(seed)

end

using .O

# %%
function mcpi_LCG(num_points = 10^9, seed = 20231226)
    seed_LCG!(seed)
    num_inside = 0
    for i in 1:num_points
        num_inside += rand_LCG()^2 + rand_LCG()^2 < 1
    end
    4num_inside / num_points
end

@time mcpi_LCG()
@time mcpi_LCG()
@time mcpi_LCG()

# %%
using LoopVectorization

isinside(i) = rand_LCG()^2 + rand_LCG()^2 < 1

function mcpi_LCG_turbo(num_points = 10^9, seed = 20231226)
    seed_LCG!(seed)
    num_inside = 0
    @turbo for i in 1:num_points
        num_inside += isinside(i)
    end
    4num_inside / num_points
end

@time mcpi_LCG_turbo()
@time mcpi_LCG_turbo()
@time mcpi_LCG_turbo()

# %%
