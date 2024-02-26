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
#     display_name: Julia 1.10.1
#     language: julia
#     name: julia-1.10
# ---

# %% tags=[]
using BenchmarkTools

function harmonic_number(n)
    s = float(zero(n))
    for k in 0:n-1
        s += 1/(n-k)
    end
    s
end

@btime harmonic_number(10^6)

# %%
harmonic_number(big(10^6))

# %%
@code_warntype harmonic_number(big(10^6))

# %%
using BenchmarkTools
using LoopVectorization

function harmonic_number_tturbo(n)
    s = float(zero(n))
    @tturbo for k in 0:n-1
        s += 1/(n-k)
    end
    s
end

@btime harmonic_number_tturbo(10^6)

# %%
using BenchmarkTools
using SpecialFunctions

harmonic_number_digamma(n) = digamma(n + 1) + MathConstants.γ

@btime harmonic_number_digamma(10^6)

# %%
@btime harmonic_number_tturbo(10^8)

# %%
@btime harmonic_number_digamma(10^8)

# %%
