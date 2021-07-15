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
a = zeros(3)
x, y, z = ntuple(i -> similar(a), Val(3))

# %%
x, y, z

# %%
x[1] = 0
x, y, z

# %%
?ntuple

# %%
