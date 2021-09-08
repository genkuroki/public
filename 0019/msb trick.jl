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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://twitter.com/mkashi/status/1434847564256268291

# %%
ENV["LINES"] = 1000

# %%
function msb64(x)
    ε = eps(Float64)
    q = (1/ε + 1) * x
    m = q - (1 - ε/2) * q
    Int64(m)
end

@show eps(Float64) == 1/2^52

[(x = 2^k; (k, log2(msb64(x - 2)), log2(msb64(x - 1)), log2(msb64(x)))) for k in 1:54]

# %%
using Quadmath

function msb128(x)
    ε = eps(Float128)
    q = (1/ε + 1) * x
    m = q - (1 - ε/2) * q
    Int128(m)
end

@show eps(Float128) == 1/Int128(2)^112

[(x = Int128(2)^k; (k, log2(msb128(x - 2)), log2(msb128(x - 1)), log2(msb128(x)))) for k in 1:114]

# %%
function msb(x)
    ε = eps(float(typeof(x)))
    q = (1/ε + 1) * x
    m = q - (1 - ε/2) * q
    oftype(x, m)
end

# %%
@show precision(Float64)
@show eps(Float64) == 1/2^(precision(Float64)-1)
[(x = 2^k; (k, log2(msb(x - 2)), log2(msb(x - 1)), log2(msb(x)))) for k in 1:precision(Float64)+1]

# %%
@show precision(BigFloat)
@show eps(BigFloat) == 1/big(2)^(precision(BigFloat)-1)
[(x = big(2)^k; (k, log2(msb(x - 2)), log2(msb(x - 1)), log2(msb(x)))) for k in 1:precision(BigFloat)+1]

# %%
