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
a = rand(1:10^6, 10^6)
length(a)

# %%
u = unique(a)
length(u)

# %%
s = sort(u)
first(s, 100) |> print

# %%
X = Set(a)
length(X)

# %%
using BenchmarkTools
R = rand(1:10^6, 100)

@btime in($a).($R)
@btime in($u).($R)
@btime Base.Fix2(insorted, $s).($R)
@btime in($X).($R);

# %%
