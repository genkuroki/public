# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.10.0-DEV
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# https://twitter.com/rtafds/status/1665041067048329219

# %%
@time using BoundaryValueDiffEq
@time using Plots
default(fmt=:png)

xspan = (0, 1)
initguess = [0.0, 0.0]
param = SciMLBase.NullParameters()

function diffeq!(dy, y, param, x)
    dy[1] = y[2] # dy₁/dx = y₂
    dy[2] = -y[1] - x # dy₂/dx = -y - x
end

function boundarycond!(residual, y, param, x)
    residual[1] = y[begin][1] # y(x_begin) = 0
    residual[2] = y[end][1] # y(x_end) = 0
end

@time bvprob = BVProblem(diffeq!, boundarycond!, initguess, xspan, param)
@time sol = solve(bvprob, GeneralMIRK4(), dt = 0.05)
@time scatter(sol.t, x -> sol(x; idxs=1); label="numerical")
@time plot!(x -> sin(x)/sin(1) - x; label="exact")

# %%
sol

# %%
?BVProblem

# %% [markdown]
# https://docs.sciml.ai/DiffEqDocs/stable/tutorials/bvp_example/

# %%
@time using BoundaryValueDiffEq
@time using Plots
default(fmt=:png)

const g = 9.81
L = 1.0
tspan = (0.0, pi / 2)
function simplependulum!(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
end

# %%
function bc1!(residual, u, p, t)
    residual[1] = u[end ÷ 2][1] + pi / 2 # the solution at the middle of the time span should be -pi/2
    residual[2] = u[end][1] - pi / 2 # the solution at the end of the time span should be pi/2
end
@time bvp1 = BVProblem(simplependulum!, bc1!, [pi / 2, pi / 2], tspan)
@time sol1 = solve(bvp1, GeneralMIRK4(), dt = 0.05)
@time plot(sol1)

# %%
