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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
function f(x::Int)::Int
    x/2
end

# %%
function g(x::Int)
    y = x/2
    z = convert(Int, y)
    Core.typeassert(Int, z)
    z
end

# %%
@code_warntype f(2)

# %%
@code_warntype g(2)

# %%
f(2)

# %%
f(3)

# %%
