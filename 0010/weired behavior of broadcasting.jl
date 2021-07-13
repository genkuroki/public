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
a = [1.0, 2.0]
b = [1.0, 2.0]
c = zeros(2)

f!(a, b, c) = @. c = b - a/b^2
@time f!(a, b, c)
@time f!(a, b, c)
@time f!(a, b, c)

# %%
a = [1.0, 2.0]
b = [1.0, 2.0]
c = zeros(2)

f2!(a, b, c) = @. c = -a/b^2 + b
@time f2!(a, b, c)
@time f2!(a, b, c)
@time f2!(a, b, c)

# %%
using BenchmarkTools
@btime f!($a, $b, $c)
@btime f2!($a, $b, $c)

# %%
g!(a, b, c) = c .= .- a ./ (b .^ 2) .+ b
@btime g!($a, $b, $c)

# %%
