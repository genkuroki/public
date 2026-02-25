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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
if haskey(ENV, "COLAB_GPU")
    import Pkg
    Pkg.add("SpecialFunctions")
    Pkg.add("CairoMakie")
end

using LinearAlgebra: norm
import CairoMakie as cm
import CairoMakie.:(..)
using SpecialFunctions
fact(x) = exp(lfactorial(x))

function vf_func(complex_func)
    function f(point)
        x, y = point
        z = complex(x, y)
        w = complex_func(z)
        w = isnan(w) ? zero(w) : w
        cm.Point2f(real(w), -imag(w))
    end
    f
end

function streamplot_complex_func(complex_func, x=-1..1, y=-1..1;
        size=(500, 500), colormap=:rainbow, color=log∘norm, title="", titlesize=20)
    fig = cm.Figure(; size)
    ax = cm.Axis(fig[1,1]; title, titlesize)
    cm.streamplot!(ax, vf_func(complex_func), x, y; colormap, color)
    fig
end

# %%
F(z) = sum(k -> z^k / fact(k), 0:5)
streamplot_complex_func(F, -2pi..2pi, -2pi..2pi)

# %%
F(z) = sum(k -> z^k / fact(k), 0:10)
streamplot_complex_func(F, -2pi..2pi, -2pi..2pi)

# %%
F(z) = sum(k -> z^k / fact(k), 0:20)
streamplot_complex_func(F, -2pi..2pi, -2pi..2pi)

# %%
F(z) = exp(z)
streamplot_complex_func(F, -2pi..2pi, -2pi..2pi)

# %%
