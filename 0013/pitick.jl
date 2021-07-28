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
using Plots, LaTeXStrings
pyplot(fmt=:svg, fontfamily="sans-serif")

function pitick(start, stop, denom; mode=:text)
    a = Int(cld(start, π/denom))
    b = Int(fld(stop, π/denom))
    tick = range(a*π/denom, b*π/denom; step=π/denom)
    ticklabel = piticklabel.((a:b) .// denom, Val(mode))
    tick, ticklabel
end

function piticklabel(x::Rational, ::Val{:text})
    iszero(x) && return "0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return S * N * "π"
    S * N * "π/" * repr(d)
end

function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return L"%$S%$N\pi"
    L"%$S\frac{%$N\pi}{%$d}"
end

a, b = -2π, 2π
plot(sin, a, b; xtick=pitick(a, b, 4), label=L"y = \sin(x)", size=(720, 250))

# %%
plot(sin, a, b; xtick=pitick(a, b, 4; mode=:latex), label=L"y = \sin(x)", size=(720, 250),
    tickfontsize=10)

# %%
gr(fmt=:auto)

plot(sin, a, b; xtick=pitick(a, b, 4), label=L"y = \sin(x)", size=(720, 250),
    fontfamily="Computer Modern", legendfontsize=12) 

# %%
plot(sin, a, b; xtick=pitick(a, b, 4; mode=:latex), label=L"y = \sin(x)", size=(720, 250),
    tickfontsize=10, fontfamily="Computer Modern", legendfontsize=12, bottom_margin=3Plots.mm)

# %%
function piticklabel(x::Rational, ::Val{:gr})
    iszero(x) && return L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return L"%$S%$N\pi"
    L"%$S\dfrac{%$N\pi}{%$d}"
end

plot(sin, a, b; xtick=pitick(a, b, 4; mode=:gr), label=L"y = \sin(x)", size=(720, 250),
    tickfontsize=10, fontfamily="Computer Modern", legendfontsize=12, bottom_margin=3Plots.mm)

# %%
