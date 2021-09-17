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
VERSION

# %%
f(x) = 3.0 * x
g(x) = x * 3.0

# %%
@code_llvm debuginfo=:none f(4.0)

# %%
@code_llvm debuginfo=:none g(4.0)

# %%
@code_native debuginfo=:none f(4.0)

# %%
@code_native debuginfo=:none g(4.0)

# %%
