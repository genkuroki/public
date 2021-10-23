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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using Plots
default(size=(400, 250))

Base.:^(f::Function, n::Integer) = n == 0 ? identity : foldr(∘, (f for _ in 1:n))

# %%
f(x) = sin(x)
plot(f, 0, 2; label="") |> display

plot([n^(1/2)*(f^n)(1.0) for n in 0:10:200]; label="")
hline!([√3]; label="")

# %%
f(x) = x - x^3/10
plot(f, 0, 2; label="") |> display

plot([n^(1/2)*(f^n)(1.0) for n in 0:10:200]; label="")

# %%
f(x) = x - x^2/4
plot(f, 0, 2; label="") |> display

plot([n*(f^n)(1.0) for n in 0:25:500]; label="")

# %%
f(x) = x - x^(3/2)/2
plot(f, 0, 2; label="") |> display

plot([n^2*(f^n)(1.0) for n in 0:25:500]; label="")

# %%
