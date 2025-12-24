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
import CairoMakie as cm
import CairoMakie.:(..)

@kwdef struct FitzhughNagumo{T}
    ϵ::T = 0.1
    s::T = 0.0
    γ::T = 1.5
    β::T = 0.8
end

param = FitzhughNagumo()

function f(point, param::FitzhughNagumo)
    x, y = point
    (; ϵ, s, γ, β) = param
    cm.Point2f((x - y - x^3 + s)/ϵ, γ*x - y + β)
end

fig, ax, pl = 
cm.streamplot(point->f(point, param), -1.5..1.5, -1.5..1.5;
    colormap = :magma)
cm.streamplot(fig[1,2], point->f(point, param), -1.5..1.5, -1.5..1.5;
    color=p->cm.RGBAf(p..., 0.0, 1.0))
fig

# %%
?cm.streamplot

# %%
import CairoMakie as cm
import CairoMakie.:(..)

@kwdef struct FitzhughNagumo{T}
    ϵ::T = 0.1
    s::T = 0.0
    γ::T = 1.5
    β::T = 0.8
end

param = FitzhughNagumo()

function f(point, param::FitzhughNagumo)
    x, y = point
    (; ϵ, s, γ, β) = param
    cm.Point2f((x - y - x^3 + s)/ϵ, γ*x - y + β)
end

fig, ax, pl = 
cm.streamplot(point->f(point, param), -1.5..1.5, -1.5..1.5;
    colormap = :magma)
#cm.streamplot(fig[1,2], point->f(point, param), -1.5..1.5, -1.5..1.5;
#    color=p->cm.RGBAf(p..., 0.0, 1.0))
fig

# %%
