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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
@inline calcM(x) = 2^x

struct MyType{N,M} end

function MyType{N}() where N
    @assert N ≥ 0
    M = calcM(N)
    MyType{N,M}()
end

@code_warntype MyType{3}() # not inferred of course

# %%
MyType{3}()

# %%
function MyType{N}() where N
    #T = typeof(N)
    M = 2^N
    MyType{N, M}()
end

@show MyType{3}()

@code_warntype MyType{3}()

# %%
@code_typed MyType{3}()

# %%
using PyCall
sys = pyimport("sys")
sys.stdout = stdout
sys.stderr = stderr

# %%
for i in 1:10
    print(i, "\n"); flush(stdout)
    sleep(0.5)
end

# %%
