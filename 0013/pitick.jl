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

# %%
using Plots

function pitick(start, stop, denom)
    a = Int(cld(start, π/denom))
    b = Int(fld(stop, π/denom))
    tick = range(a*π/denom, b*π/denom; step=π/denom)
    ticklabel = piticklabel.(a:b, denom)
    tick, ticklabel
end

function piticklabel(x, denom)
    d = repr(denom)
    x == 0 && return "0"
    x == 1 && return "π/" * d
    x == -1 && return "-π/" * d
    if mod(x, denom) == 0
        q = x ÷ denom
        q == 1 && return "π"
        q == -1 && return "-π"
        return repr(q) * "π"
    end
    return repr(x) * "π/" * d
end

a, b = -2π, 2π
plot(sin, a, b; xtick=pitick(a, b, 4), label="y = sin(x)", size=(720, 300))

# %%
