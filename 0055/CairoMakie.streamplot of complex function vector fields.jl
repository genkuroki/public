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
streamplot_complex_func(z->one(z), -1..1, -1..1; title="f(z) = 1")

# %%
streamplot_complex_func(z->z, -1..1, -1..1; title="f(z) = z")

# %%
streamplot_complex_func(z->z^2, -1..1, -1..1; title="f(z) = z²")

# %%
streamplot_complex_func(z->z^3, -1..1, -1..1; title="f(z) = z³")

# %%
streamplot_complex_func(z->z^4, -1..1, -1..1; title="f(z) = z⁴")

# %%
streamplot_complex_func(z->z*(z-1), -2..3, -2.5..2.5; title="f(z) = z(z-1)")

# %%
streamplot_complex_func(z->1/z, -1..1, -1..1; title="f(z) = 1/z")

# %%
streamplot_complex_func(z->1/z^2, -1..1, -1..1; title="f(z) = 1/z²")

# %%
streamplot_complex_func(z->1/z^3, -1..1, -1..1; title="f(z) = 1/z³")

# %%
streamplot_complex_func(z->1/z^4, -1..1, -1..1; title="f(z) = 1/z⁴")

# %%
streamplot_complex_func(z->1/(z+0.5)-1/(z-0.5), -1..1, -1..1; 
    title="f(z) = 1/(z+0.5) - 1/(z-0.5)")

# %%
streamplot_complex_func(z->1/(z+0.1) - 1/(z-0.1), -1..1, -1..1; 
    title="f(z) = 1/(z+0.1) - 1/(z-0.1)")

# %%
streamplot_complex_func(z->2/(z+0.5) - 1/z - 1/(z-0.5), -1..1, -1..1;
    title="f(z) = 2/(z+0.5) - 1/z - 1/(z-0.5)")

# %%
streamplot_complex_func(z->2/(z+0.1) - 1/z - 1/(z-0.1), -1..1, -1..1;
    title="f(z) = 2/(z+0.1) - 1/z - 1/(z-0.1)")

# %%
streamplot_complex_func(z->1/(im*z), -2..2, -2..2; title="f(z) = 1/(iz)")

# %%
streamplot_complex_func(z->(1+2im)/z, -2..2, -2..2; title="f(z) = (1+2i)/z")

# %%
streamplot_complex_func(z->(1+2im)/(z+0.5) - (1+2im)/(z-0.5), -2..2, -2..2;
    title="f(z) = (1+2i)/(z+0.5) - (1+2i)/(z-0.5)")

# %%
f(z) = (1+2im)/(z+1+im) - (1+2im)/(z-1-im) + (1+2im)/(z+1-im) - (1+2im)/(z-1+im)
title = "f(z) = (1+2i)/(z+1+i) - (1+2i)/(z-1-i)\n       + (1+2i)/(z+1-i) - (1+2i)/(z-1+i)"
streamplot_complex_func(f, -2.5..2.5, -2.5..2.5, size=(500, 520); title)

# %%
streamplot_complex_func(exp, -2..2, -2..2; title="f(z) = exp(z)")

# %%
streamplot_complex_func(z->exp(1/z), -0.3..0.3, -0.3..0.3; title="f(z) = exp(1/z)")

# %%
streamplot_complex_func(sinpi, -2..2, -2..2; title="f(z) = sin(πz)")

# %%
streamplot_complex_func(z->1/sinpi(z), -2..2, -2..2; title="f(z) = 1/sin(πz)")

# %%
streamplot_complex_func(tanpi, -2..2, -1..1; title="f(z) = tan(πz)")

# %%
streamplot_complex_func(gamma, -8..8, -5..5; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(gamma, -3..1, -1..1; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(gamma, -7..1, -2..2; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(gamma, 0..3, -1..1; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(gamma, 0..8, -2..2; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(gamma, 0..16, -4..4; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(gamma, 1..3, 0..8; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(gamma, 1..5, 0..16; title="f(z) = Γ(z)")

# %%
streamplot_complex_func(zeta, 0..2, -1..1; title="f(z) = ζ(z) at z = 1", size=(400, 400))

# %%
streamplot_complex_func(zeta, -5..5, -5..5; title="f(z) = ζ(z)", size=(400, 400))

# %%
streamplot_complex_func(zeta, -5..5, 10..20; title="f(z) = ζ(z)", size=(400, 400))

# %%
streamplot_complex_func(zeta, 0.2..0.8, 13..15; 
    title="f(z) = ζ(z) at the first nontrivial zero", size=(400, 400))

# %%
streamplot_complex_func(zeta, 0.2..0.8, 20..22; 
    title="f(z) = ζ(z) at the second nontrivial zero", size=(400, 400))

# %%
streamplot_complex_func(zeta, 0.2..0.8, 24..26; 
    title="f(z) = ζ(z) at the third nontrivial zero", size=(400, 400))

# %%
