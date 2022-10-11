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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
Int

# %%
promote_type(Int, Float64)

# %%
f(x) = 2x
@code_typed f(3.0)

# %%
promote_type(Int, Float32)

# %%
@code_typed f(3f0)

# %%
promote_type(Int, ComplexF64)

# %%
@code_typed f(complex(3.0))

# %%
promote_type(Int, Rational{Int})

# %%
@code_typed f(3//1)

# %%
