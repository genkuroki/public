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
using Roots

function withcounter(f, n = 0)
    c = Ref(n)
    F(x...; y...) = (c[] += 1; f(x...; y...))
    getcounter() = c[]
    setcounter!(n) = c[] = n
    F, getcounter, setcounter!
end

f(x::Float64)::Float64 = @ccall sin(x::Float64)::Float64
@show f(1.0);

# %%
F, getcounter, setcounter! = withcounter(x -> f(x) - 0.8414709848078965)
@show find_zero(F, 0.9)
@show getcounter();

# %%
setcounter!(0)
@show find_zeros(F, 0, 2π)
@show getcounter();

# %%
dump(F)

# %%
