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

export TaskLocalXorshift64, XOS

using Random: Random, AbstractRNG, RandomDevice, SamplerType, SamplerTrivial, CloseOpen01_64

struct TaskLocalXorshift64 <: AbstractRNG end

const XOS = TaskLocalXorshift64()

@inline getstate(::TaskLocalXorshift64) = current_task().rngState0

@inline function setstate!(xos::TaskLocalXorshift64, seed::Integer)
    current_task().rngState0 = mod(seed, UInt64)
    xos
end

@inline function Random.rand(xos::TaskLocalXorshift64, ::SamplerType{UInt64})
    x = res = getstate(xos)
    x ⊻= x << 13
    x ⊻= x >> 7
    x ⊻= x << 17
    setstate!(xos, x)
    res
end

@inline function Random.rand(xos::TaskLocalXorshift64, ::SamplerTrivial{CloseOpen01_64})
    (rand(xos, UInt64) & UInt64(2^52 - 1)) / 2^52
end

Random.seed!(xos::TaskLocalXorshift64) =
    setstate!(xos, rand(RandomDevice(), UInt64))

Random.seed!(xos::TaskLocalXorshift64, seed::Integer) =
    setstate!(xos, seed)

end

using .O

@show O.getstate(XOS)
A = [rand(XOS) for _ in 1:20]

# %%
Random.seed!(XOS)
@show O.getstate(XOS)
Random.seed!(XOS, 0x1234567)
@show O.getstate(XOS);

# %%
function mcpi(L=10^9, rng=XOS)
    c = 0
    for i in 1:L
        c += rand(rng)^2 + rand(rng)^2 ≤ 1
    end
    4c/L
end

@time mcpi()
@time mcpi()
@time mcpi()

# %%
using LoopVectorization

@inline isinside(i, rng=XOS) =
    rand(rng)^2 + rand(rng)^2 ≤ 1

function mcpi_turbo(L=10^9, isinsiderng=isinside)
    c = 0
    @turbo for i in 1:L
        c += isinsiderng(i)
    end
    4c/L
end

@time mcpi_turbo()
@time mcpi_turbo()
@time mcpi_turbo()

# %%
function mcpi_tturbo(L=10^9, isinsiderng=isinside)
    c = 0
    @tturbo for i in 1:L
        c += isinsiderng(i)
    end
    4c/L
end

mcpi_tturbo()
@time mcpi_tturbo()
@time mcpi_tturbo()
@time mcpi_tturbo()

# %%
