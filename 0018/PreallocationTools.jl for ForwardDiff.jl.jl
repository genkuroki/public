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
# https://discourse.julialang.org/t/how-to-make-this-function-compatible-with-forwarddiff/67415

# %%
using ForwardDiff
using PreallocationTools

function f(x, tmp)
    tmp = get_tmp(tmp, x)
    @. tmp = x + 1
    sum(tmp)
end

x = [1.0, 2.0, 3.0]
tmp = dualcache(similar(x))
@show f(x, tmp)
ForwardDiff.gradient(x -> f(x, tmp), x)

# %%
