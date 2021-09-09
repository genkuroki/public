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
using LinearAlgebra

S(x, A) = dot(x, A, x)
S(x) = S(x, A)

niters = 10^6
n = 10
A = randn(n, n)
x = randn(n);

# %%
@code_warntype S(x, A)

# %%
@code_warntype S(x)

# %%
function f(niters, x, A)
    y = Vector{eltype(x)}(undef, niters)
    for i in 1:niters
        y[i] = S(x, A)
    end
    y
end

@code_warntype f(niters, x, A)

# %%
function f(niters, x)
    y = Vector{eltype(x)}(undef, niters)
    for i in 1:niters
        y[i] = S(x)
    end
    y
end

@code_warntype f(niters, x)

# %%
using BenchmarkTools

@btime f(niters, $x, $A)
@btime f(niters, $x);

# %%
