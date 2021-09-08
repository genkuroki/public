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

function msb64(x)
    ε = eps(Float64)
    q = (1/ε + 1) * x
    m = q - (1 - ε/2) * q
    Int64(m)
end

[(x = 2^k; (k, log2(msb64(x - 2)), log2(msb64(x - 1)), log2(msb64(x)))) for k in 1:54]

# %%
using Quadmath

function msb128(x)
    ε = eps(Float128)
    q = (1/ε + 1) * x
    m = q - (1 - ε/2) * q
    Int128(m)
end

[(x = Int128(2)^k; (k, log2(msb128(x - 2)), log2(msb128(x - 1)), log2(msb128(x)))) for k in 1:114]

# %%
