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
f(x::Float64) = 2.0x + 1.0
@code_llvm debuginfo=:none f(1.2)

# %%
g(x) = 2x + 1
@code_llvm debuginfo=:none g(1.2)

# %%
?@code_llvm

# %%
