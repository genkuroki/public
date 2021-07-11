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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using LinearAlgebra, BenchmarkTools

A = randn(10, 1000)
B = randn(1000, 1000)
C = randn(1000, 10000)
P = zeros(10, 1000)
Q = zeros(10, 10000)
R = zeros(1000, 10000)
S = zeros(10, 10000)

"""(A * B) * C"""
function f!(A, B, C, P, Q)
    mul!(P, A, B)
    mul!(Q, P, C)
end

"""A * (B * C)"""
function g!(A, B, C, R, S)
    mul!(R, B, C)
    mul!(S, A, R)
end

f!(A, B, C, P, Q) ≈ g!(A, B, C, R, S)

# %%
@benchmark f!($A, $B, $C, $P, $Q)

# %%
@benchmark g!($A, $B, $C, $R, $S)

# %%
