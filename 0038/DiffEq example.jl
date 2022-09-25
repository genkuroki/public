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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
using DifferentialEquations
using PyPlot: PyPlot, plt

# %%
prob = ODEProblem((u, p, t) -> u*cos(t), 1.0, (0.0, 2.0))
sol = solve(prob)

# %%
f(t) = exp(sin(t))
t = range(sol.prob.tspan..., 1000)
plt.figure(figsize=(6, 4))
plt.plot(t, f.(t); label="analytical solution")
plt.plot(t, sol.(t); label="numerical solution", ls="--")
plt.legend()
plt.xlabel("t")
plt.ylabel("u(t)")

# %%
