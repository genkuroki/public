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

# %% [markdown]
# [Run on Google Colab](https://colab.research.google.com/github/genkuroki/public/blob/main/0055/CairoMakie.streamplot%20of%20complex%20function%20vector%20fields.ipynb)

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

function vf_func(complex_func)
    function f(point)
        x, y = point
        z = complex(x, y)
        cm.Point2f(real(complex_func(z)), -imag(complex_func(z)))
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

function streamplot_complex_func(complex_func_str::AbstractString,  x=-1..1, y=-1..1;
        size=(500, 500), colormap=:rainbow, color=log∘norm, title="", titlesize=20)
    title == "" && (title = "f(z) = " * complex_func_str)
    @eval f(z) = $(Meta.parse(complex_func_str))
    streamplot_complex_func(p->invokelatest(f, p), x, y; size, colormap, color, title, titlesize)
end

# %%
streamplot_complex_func("1", -1..1, -1..1)

# %%
streamplot_complex_func("z", -1..1, -1..1)

# %%
streamplot_complex_func("z^2", -1..1, -1..1)

# %%
streamplot_complex_func("z^3", -1..1, -1..1)

# %%
streamplot_complex_func("z^3", -1..1, -1..1)

# %%
streamplot_complex_func("1/z", -1..1, -1..1)

# %%
streamplot_complex_func("1/z^2", -1..1, -1..1)

# %%
streamplot_complex_func("1/z^3", -1..1, -1..1)

# %%
streamplot_complex_func("1/(z+0.5) - 1/(z-0.5)", -1..1, -1..1)

# %%
streamplot_complex_func("1/(z+0.1) - 1/(z-0.1)", -1..1, -1..1)

# %%
streamplot_complex_func("2/(z+0.5) - 1/z - 1/(z-0.5)", -1..1, -1..1)

# %%
streamplot_complex_func("2/(z+0.1) - 1/z - 1/(z-0.1)", -1..1, -1..1)

# %%
streamplot_complex_func("im/z", -2..2, -2..2)

# %%
streamplot_complex_func("(1+2im)/z", -2..2, -2..2)

# %%
streamplot_complex_func("(1+2im)/(z+0.5) - (1+2im)/(z-0.5)", -2..2, -2..2)

# %%
complex_func_str = "(1+2im)/(z+1+im) - (1+2im)/(z-1-im) +\n (1+2im)/(z+1-im) - (1+2im)/(z-1+im)"
streamplot_complex_func(complex_func_str, -2.5..2.5, -2.5..2.5, size=(500, 520))

# %%
