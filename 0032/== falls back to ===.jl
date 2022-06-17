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
module O
struct Foo{T} a::T end
end

# %%
foo = O.Foo(1)
bar = O.Foo(2)

# %%
@which foo == bar

# %% [markdown]
# ```julia
# ==(x, y) = x === y
# ```
#
# See https://docs.julialang.org/en/v1/base/math/#Base.:==

# %%
