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
@show all(f, 1:10^8)

g(k) = highestdigit(big(2)^k) == highestdigit_log10(big(2)^k)
@show all(g, 1:10^4)

[(n, highestdigit(n)) for n in (rand(1:10^8) for _ in 1:20)]

# %%
highestdigit_naive(x::Integer) = x ÷ 10^(ndigits(x) - 1)

@show x = Int128(2)^102
@show highestdigit_naive(x)
@show highestdigit(x);

# %%
function highestdigit_log10_naive1(x)
    x == 0 && return 0
    r = log10(x) % 1
    for k in 1:8
        r < log10(oftype(x, k+1)) && return k
    end
    9
end

@show x = 30
@show highestdigit_log10_naive1(x)
@show highestdigit_log10(x)
@show e = eps(float(one(x)))
@show log10(x) % 1
@show log10(x) % 1 + e
@show log10(oftype(x, 3))
@show log10(x) % 1 < log10(oftype(x, 3))
@show log10(x) % 1 + e < log10(oftype(x, 3));

# %%
function highestdigit_log10_naive2(x; e=eps(float(one(x))))
    x == 0 && return 0
    r = log10(x) % 1
    for k in 1:8
        r + e < log10(k+1) && return k
    end
    9
end

@show x = big(2)
@show highestdigit_log10_naive2(x)
@show highestdigit_log10(x)
@show e = eps(float(one(x)))
@show log10(x) % 1 + e
@show log10(1 + 1)
@show log10(oftype(x, 1 + 1))
@show log10(x) % 1 + e < log10(1 + 1)
@show log10(x) % 1 + e < log10(oftype(x, 1 + 1));

# %%
