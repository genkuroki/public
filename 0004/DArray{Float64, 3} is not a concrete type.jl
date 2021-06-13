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
#     display_name: Julia 1.7.0-DEV
#     language: julia
#     name: julia-1.7
# ---

# %%
using BenchmarkTools
using DistributedArrays
struct Foo{T, N} d::DArray{T, N} end
struct Bar{T, N, A} d::DArray{T, N, A} end

# %%
d = drand(4, 3, 2)
l = localpart(d)
foo = Foo(d)
bar = Bar(d);

# %%
f(foobar, i, j, k) = localpart(foobar.d)[i, j, k]
g(d, i, j, k) = localpart(d)[i, j, k]

# %%
@btime f($foo, 3, 2, 1)

# %%
@btime f($bar, 3, 2, 1)

# %%
@btime g($d, 3, 2, 1)

# %%
@btime getindex($l, 3, 2, 1)

# %%
@code_warntype f(foo, 3, 2, 1)

# %%
@code_warntype f(bar, 3, 2, 1)

# %%
isconcretetype(Array{Float64, 3})

# %%
isconcretetype(DArray{Float64, 3})

# %%
typeof(d)

# %%
isconcretetype(DArray{Float64, 3, Array{Float64, 3}})

# %%
