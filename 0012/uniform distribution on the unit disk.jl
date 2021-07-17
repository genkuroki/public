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
using Plots

function randunitdisk()
    y, x = sincospi(2rand())
    r = √rand()
    r*x, r*y
end

randunitdisk(n) = [randunitdisk() for _ in 1:n]

# %%
scatter(randunitdisk(2^12); legend=false, msw=0, ms=3, alpha=0.5, size=(500, 500))

# %%
scatter(randunitdisk(2^14); legend=false, msw=0, ms=1.5, alpha=0.5, size=(500, 500))

# %%
