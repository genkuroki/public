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
methods(factorial)

# %%
factorial(20)

# %%
factorial(21)

# %%
@which factorial(21)

# %%
Base._fact_table64

# %%
Base._fact_table64[1] = 666

# %%
factorial(1)

# %%
factorial(0.5)

# %%
using SpecialFunctions

# %%
factorial(0.5)

# %%
gamma(1.5)

# %%
@which factorial(0.5)

# %%
methods(factorial)

# %%
factorial(20)

# %%
@which factorial(20)

# %%
