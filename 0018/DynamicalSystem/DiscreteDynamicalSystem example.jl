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
using DynamicalSystems

# %%
f(u, p, t) = SVector(u[2], u[1]+u[2])
u0 = [0, 1]
prob = DiscreteDynamicalSystem(f, u0, nothing)

# %% tags=[]
tr = trajectory(prob, 20, Ttr=10)

# %% tags=[]
tr = trajectory(prob, 20, Ttr=0)

# %%
dump(tr)

# %%
