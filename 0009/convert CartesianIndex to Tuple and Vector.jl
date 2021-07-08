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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %%
ci = CartesianIndex(2, 3, 4)

# %%
dump(ci)

# %%
Tuple(ci)

# %%
collect(Tuple(ci))

# %%
