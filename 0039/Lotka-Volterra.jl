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
#     display_name: Julia 1.8.4
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://gendai.media/articles/-/63904?page=3
#
# ![2022-12-30.png](attachment:63251af2-2ec0-43af-92a0-5cf04ab6868f.png)
#
# https://twitter.com/kimu3_slime/status/1608678692817080320?s=20&t=DN3hlMCpGLP6Q9HeHs8gaQ

# %%
using OrdinaryDiffEq
using StaticArrays
using Plots
default(fmt=:png)

function LotkaVolterra(u, param, t)
    S, W = u
    p, q, r, s = param
    dS =  p*S - q*W*S
    dW = -r*W + s*W*S
    SVector(dS, dW)
end

param = [1.4, 0.7, 1.0, 1.0]
u0 = SVector(2.0, 1.0)
tspan = (0, 25)
prob = ODEProblem(LotkaVolterra, u0, tspan, param)

# %%
sol = solve(prob, Euler(); dt=0.05)
plot(sol; legend=false, ls=[:solid :dash])

# %%
sol = solve(prob, Vern7())
plot(sol; legend=false, ls=[:solid :dash])

# %%
