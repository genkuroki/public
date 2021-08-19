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
A = 10(1:9) .+ (1:9)'

# %%
i = 4
@views A00 = A[begin:i-1, begin:i-1]

# %%
@views a10 = A[begin:i-1, i:i]

# %%
@views A02 = A[begin:i-1, i+1:end]

# %%
@views A10 = A[i:i, begin:i-1]

# %%
@views a11 = A[i:i, i:i]

# %%
@views a12 = A[i:i, i+1:end]

# %%
@views A20 = A[i+1:end, begin:i-1]

# %%
@views a21 = A[i+1:end, i:i]

# %%
@views A22 = A[i+1:end, i+1:end]

# %%
