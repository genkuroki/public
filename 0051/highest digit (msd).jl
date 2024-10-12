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
#     display_name: Julia 1.11.0
#     language: julia
#     name: julia-1.11
# ---

# %%
highestdigit(x::Integer) = x ÷ oftype(x, 10)^(ndigits(x) - 1)

function highestdigit_log10(x; e=eps(float(one(x))))
    x == 0 && return 0
    r = log10(x) % 1
    for k in 1:8
        r + e < log10(oftype(x, k+1)) && return k
    end
    9
end

f(k) = highestdigit(k) == highestdigit_log10(k)
g(k) = highestdigit(big(2)^k) == highestdigit_log10(big(2)^k)
@show all(f, 1:10^8)
@show all(g, 1:10^4)

[(n, highestdigit(n)) for n in (rand(1:10^8) for _ in 1:20)]

# %%
