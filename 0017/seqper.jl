# -*- coding: utf-8 -*-
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

# %% [markdown]
# https://discourse.julialang.org/t/julia-implementation-of-seqperiod-from-matlab/66791

# %%
using BenchmarkTools

function testvec(n, per; tol=1e-3)
    X = [float(mod1(k, per)) for k in 1:n]
    @. X + tol * (rand() - 0.5)
end

X = testvec(10^4, 999)
a = zeros(length(X))
a[end] = 2e-3
Y = X + a;

# %%
function seqper(x; tol=0.001)
    ind1 = findall(≤(tol), [abs.(x .- x[1])...])
    period = length(x)
    for i = 2:length(ind1)
        if maximum(abs.(x[ind1[i]:end] .- x[1:end-ind1[i]+1])) ≤ tol
            period = ind1[i] - 1
            break
        end
    end 
    return period
end

@show seqper(X)
@btime seqper($X)
@show seqper(Y)
@btime seqper($Y);

# %%
function seqper1(x; tol=1e-3)
    @inbounds for k in 2:length(x)
        if abs(x[k] - x[1]) ≤ tol
            all(j -> abs(x[j] - x[j-k+1]) ≤ tol, k:lastindex(x)) && return k - 1
        end
    end 
    return length(x)
end

@show seqper1(X)
@btime seqper1($X)
@show seqper1(Y)
@btime seqper1($Y);

# %%
