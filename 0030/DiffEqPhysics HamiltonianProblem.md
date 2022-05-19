---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Julia 1.7.2
    language: julia
    name: julia-1.7
---

https://discourse.julialang.org/t/how-to-define-and-integrate-a-hamiltonianproblem/81303

<!-- #region -->
```julia
using DiffEqPhysics
using DifferentialEquations

function cart_hamiltonian(q, p; g = 9.8, m1 = 2.0 , m2 = 1.0, L = 0.5) 
    (x, θ) = q
    (p_x, p_θ) = p
    h =  g * L * m2 * cos(θ) + 
        (L^2 * m2 * p_x^2 + (m1 + m2)* p_θ^2 - 
        2 *L * m2 * p_x * p_θ * cos(θ) ) / 
        ( 2 * L^2*m2*(m1 + m2 * sin(θ)^2))
    return h
end

a = HamiltonianProblem(cart_hamiltonian,[0,0.3],[0,.1],[0,.2]) # gives me ODE problem?
sol = init(a, Tsit5())
```

```
┌ Warning: Hamiltonians with 3 arguments are deprecated; please use `H(p, q, params, t)`
│   caller = #HamiltonianProblem#1 at hamiltonian.jl:30 [inlined]
└ @ Core D:\.julia\packages\DiffEqPhysics\SeVwO\src\hamiltonian.jl:30
MethodError: no method matching cart_hamiltonian(::Vector{Float64}, ::Vector{ForwardDiff.Dual{ForwardDiff.Tag{DiffEqPhysics.PhysicsTag, Float64}, Float64, 2}}, ::Nothing)
Closest candidates are:
  cart_hamiltonian(::Any, ::Any; g, m1, m2, L) at In[1]:4
```
<!-- #endregion -->

```julia
module Revision1

using DiffEqPhysics
using DifferentialEquations
using Plots

function cart_hamiltonian(p, q, params, t; g = 9.8, m1 = 2.0 , m2 = 1.0, L = 0.5) 
    (x, θ) = q
    (p_x, p_θ) = p
    h =  g * L * m2 * cos(θ) + 
        (L^2 * m2 * p_x^2 + (m1 + m2)* p_θ^2 - 
        2 *L * m2 * p_x * p_θ * cos(θ) ) / 
        ( 2 * L^2*m2*(m1 + m2 * sin(θ)^2))
    return h
end

a = HamiltonianProblem(cart_hamiltonian,[0,0.3],[0,.1],[0,.2]) # gives me ODE problem?
sol = solve(a, Tsit5())
plot(sol; legend=:outertopright) |> display

end
```

```julia
module Revision2

using DiffEqPhysics
using DifferentialEquations
using Plots
using StaticArrays

params = (g = 9.8, m1 = 2.0 , m2 = 1.0, L = 0.5)

function cart_hamiltonian(p, q, params, t)
    (; g, m1, m2, L) = params
    (x, θ) = q
    (p_x, p_θ) = p
    h =  g * L * m2 * cos(θ) + 
        (L^2 * m2 * p_x^2 + (m1 + m2)* p_θ^2 - 
        2 *L * m2 * p_x * p_θ * cos(θ) ) / 
        ( 2 * L^2*m2*(m1 + m2 * sin(θ)^2))
    return h
end

p0 = SVector(0.0, 0.3)
q0 = SVector(0.0, 0.1)
tspan = (0.0, 0.2)
a = HamiltonianProblem(cart_hamiltonian, p0, q0, tspan, params)
sol = solve(a, Tsit5())
plot(sol; legend=:outertopright) |> display

end
```

```julia
Revision2.sol
```

```julia

```
