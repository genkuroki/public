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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Nemo

# %%
F4, ω = finite_field(4, "ω")

# %%
ω^2

# %%
F4x, x = F4["x"]

# %%
f = x^3 + x^2 + x + ω

# %%
F64, α = finite_field(f, "α")-

# %%
order(F64)

# %%
ENV["LINES"] = 200
[(k, α^k) for k in 1:63]

# %%
(x^3 + x^2 + x + ω) * (x^3 + x^2 + x + ω + 1)

# %%
F2, _ = finite_field(2)

# %%
F2t, t = F2["t"]

# %%
g = t^6 + t^4 + t^3 + t + 1

# %%
FF64, a = finite_field(g, "a")

# %%
ENV["LINES"] = 200
[(k, a^k) for k in 1:63]

# %%
