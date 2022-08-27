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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://twitter.com/dannchu/status/1563401132902588417
#
# ![FbJQ9oxacAAPvN2.png](attachment:120d12de-b14e-4749-a52f-09160a479bd4.png)

# %%
using Optim

function F(x; a=0.39, b=0.30, c=0.21, d=0.10)
    p, q, r = x
    (p^2 + 2p*r - a)^2 + (r^2 - b)^2 + (q^2 + 2q*r - c)^2 + (2p*q - d)^2
end

o = optimize(F, fill(1/3, 3), LBFGS())

# %%
o.minimizer

# %%
sum(o.minimizer)

# %%
