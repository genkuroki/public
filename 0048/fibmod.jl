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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
function fibmod(k, N)
    F = zeros(Int, k)
    F[1] = F[2] = 1
    for i in 1:k-2
        F[i+2] = mod(F[i+1] + F[i], N)
    end
    F
end

F = fibmod(2024, 13)
@show F[2024]
@show F[1:28]
@show F[29:56];

# %%
using Primes

# \boxplus TAB → ⊞
a ⊞ b = length(primes(prime(a) + prime(b)))

for (a, b) in ((1, 1), (1, 2), (3, 4), (5, 2), (5, 7), (6, 9), (10, 20), (100, 200))
    @eval @show $a ⊞ $b
end

# %%
