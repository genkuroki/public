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
(((x, y),) -> x + y).([(1, 2), (3, 4)])

# %%
f((x, y)) = x + y
f.([(1, 2), (3, 4)])

# %%
(t -> +(t...)).([(1, 2), (3, 4)])

# %%
