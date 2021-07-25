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
# https://discourse.julialang.org/t/problem-with-distrubuted-package/65199/5

# %%
using Distributed
using DifferentialEquations
# using Plots

@show addprocs()

@everywhere using DifferentialEquations
# Linear ODE which starts at 0.5 and solves from t=0.0 to t=1.0
@everywhere f(u,p,t) = 1.01*u
prob = ODEProblem(f,0.5,(0.0,1.0))

@everywhere function prob_func(prob,i,repeat) # i is the unique id 1:trajectories for the problem
  remake(prob,u0=rand()*prob.u0)
end

ensemble_prob = EnsembleProblem(prob,prob_func=prob_func)
sim = solve(ensemble_prob,Tsit5(),EnsembleDistributed(),trajectories=100)

# %%
versioninfo()

# %%
using Pkg
Pkg.status("DifferentialEquations")
Pkg.status("DiffEqBase"; mode=PKGMODE_MANIFEST)
Pkg.status("SciMLBase"; mode=PKGMODE_MANIFEST)

# %%
