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
f(x::Int)::Int = x > 0 ? x / 2 : x
g(x::Int) = x > 0 ? typeassert(convert(Int, x / 2), Int) : typeassert(convert(Int, x), Int)

# %%
@code_warntype f(1)

# %%
@code_warntype g(1)

# %%
@code_typed f(1)

# %%
@code_typed g(1)

# %%
