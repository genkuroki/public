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
using Random

module O

export TaskLocalLCG, LCG

using Random: Random, AbstractRNG, RandomDevice

const a = Int32(48271)
const m = Int32(2^31 - 1)
const q = m ÷ a
const r = m % a

struct TaskLocalLCG <: AbstractRNG end
const LCG = TaskLocalLCG()

@inline getstate(::TaskLocalLCG) = mod(current_task().rngState0, Int32)
@inline setstate!(lcg::TaskLocalLCG, x::Integer) = (current_task().rngState0 = x; lcg)

@inline function Random.rand(lcg::TaskLocalLCG) # 手抜き
    x = getstate(lcg)
    hi, lo = divrem(x, q)
    x = a * lo - r * hi
    x = ifelse(x > zero(x), x, x + m)
    setstate!(lcg, x)
    x / m # 本当は (x - 1) / m とするべき。
end

Random.seed!(lcg::TaskLocalLCG) = setstate!(lcg, rand(RandomDevice(), UInt64))
Random.seed!(lcg::TaskLocalLCG, x) = setstate!(lcg, x)

end

using .O

# %%
function mcpi_LCG(num_points = 10^9, seed = 20231226)
    Random.seed!(LCG, seed)
    num_inside = 0
    for i in 1:num_points
        num_inside += rand(LCG)^2 + rand(LCG)^2 < 1
    end
    4num_inside / num_points
end

@time mcpi_LCG()
@time mcpi_LCG()
@time mcpi_LCG()
@time mcpi_LCG()
@time mcpi_LCG()
@time mcpi_LCG()

# %%
using LoopVectorization

@inline isinside(i) = rand(LCG)^2 + rand(LCG)^2 < 1

function mcpi_LCG_turbo(num_points = 10^9, seed = 20231226)
    Random.seed!(LCG, seed)
    num_inside = 0
    @turbo for i in 1:num_points
        num_inside += isinside(i)
    end
    4num_inside / num_points
end

@time mcpi_LCG_turbo()
@time mcpi_LCG_turbo()
@time mcpi_LCG_turbo()
@time mcpi_LCG_turbo()
@time mcpi_LCG_turbo()
@time mcpi_LCG_turbo()

# %%
function mcpi_LCG_tturbo(num_points = 10^9, seed = 20231226)
    Random.seed!(LCG, seed)
    num_inside = 0
    @tturbo for i in 1:num_points
        num_inside += isinside(i)
    end
    4num_inside / num_points
end

mcpi_LCG_tturbo()
@time mcpi_LCG_tturbo()
@time mcpi_LCG_tturbo()
@time mcpi_LCG_tturbo()
@time mcpi_LCG_tturbo()
@time mcpi_LCG_tturbo()
@time mcpi_LCG_tturbo()

# %%
