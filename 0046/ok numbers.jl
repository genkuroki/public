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

# %% [markdown]
# https://x.com/dannchu/status/1740282046675894744

# %%
function d(n)
    @assert n > 0
    r = one(n)
    for k in 2:isqrt(n)
        e = 0
        while n % k == 0
            e += 1
            n ÷= k
        end
        r *= e + 1
    end
    n > 1 ? 2r : r
end

issquare(n) = isqrt(n)^2 == n

function isok(n)
    n == 1 && return false
    issquare(n+1) && issquare(d(n)) && issquare(isqrt(d(n))) && issquare(d(d(n-1)))
end

function list_of_ok_numbers(L)
    [n for n in 1:L if isok(n)]
end

@time list_of_ok_numbers(10^5)
@time list_of_ok_numbers(10^5)
@time list_of_ok_numbers(10^5)

# %%
