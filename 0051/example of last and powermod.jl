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
#     display_name: Julia 1.11.1
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# https://x.com/dannchu/status/1855466229332594870

# %%
using BenchmarkTools

f() = maximum(n for n in 11:99 if powermod(3, n, 5) == 1)
g() = first(n for n in 99:-1:11 if powermod(3, n, 5) == 1)
h() = last(n for n in 11:99 if powermod(3, n, 5) == 1)

@btime f()
@btime g()
@btime h()
@show(f()) == @show(g()) == @show(h())

# %%
gen = (n for n in 11:99 if powermod(3, n, 5) == 1)

# %%
last(gen)

# %%
@which last(gen)

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/8f5b7ca12ad48c6d740e058312fc8cf2bbe67848/base/generator.jl#L56

# %%
gen.iter

# %%
last(gen.iter)

# %%
@which last(gen.iter)

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/8f5b7ca12ad48c6d740e058312fc8cf2bbe67848/base/iterators.jl#L561

# %%
Iterators.reverse(gen.iter)

# %%
first(Iterators.reverse(gen.iter))

# %%
