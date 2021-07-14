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

# %% [markdown] tags=[]
# Optimization by specialization to argument types comes with the tradeoff of increased compilation time. An extreme case of this is shown in the following example.

# %%
f(x) = x^2 + 3x + 2
g(::Val{x}) where x = f(x)
f(10^6), g(Val(10^6))

# %%
using BenchmarkTools
@btime f(10^6)
@btime g($(Val(10^6)));

# %% [markdown]
# `g(Val(10^6))` is ultra fast because it is specialized to the argument type `Val{10^6}` and compiled to `return 1000003000002`.

# %%
@code_typed debuginfo=:none g(Val(10^6))

# %% [markdown]
# However, `g(Val(k))` is compiled separately for each different `k`.  So if you run `g(Val(k))` for a large number of different `k`'s,  it will perform a large number of compilations and will be very slow. (After compilation, though, it will be explosively fast.) 

# %%
F(n) = [f(k) for k in 1:n]
G(n) = [g(Val(k)) for k in 1:n]

@time F(10^4)
@time G(10^4);

# %% [markdown]
# The first execution of `G(10^4)` is very slow.

# %%
F(10^4) == G(10^4)

# %% [markdown]
# Thus, it is not reasonable to try to optimize by specialization to the argument types by different large numbers of `Val{k}` types.
#
# On the other hand, the native code specialized to argument types is very fast, so if compilation time is not an issue, optimization by specialization to argument types should be done aggressively.
#
# In short, it is a matter of trade-off.

# %%
