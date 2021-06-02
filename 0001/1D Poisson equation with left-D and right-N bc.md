---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Julia 1.7.0-DEV
    language: julia
    name: julia-1.7
---

See

* https://twitter.com/OPPO89694572/status/1399910007777349633
* https://twitter.com/OPPO89694572/status/1399721433127923722

```julia
using LinearAlgebra
using Plots
using Printf
```

```julia
"""
1-dimensional discrete Laplacian
with left Dirichlet and right Neumann boundary conditions

Laplacian(N) × (-1) is equal to the Cartan matrix of B-type.
"""
Laplacian(N) = Tridiagonal(
    [fill(1, N-2); 2],
    [fill(-2, N-1); -2],
    fill(1, N-1)
)

BC(N, h, α, β) = [α; zeros(N-2); 2β*h]

Laplacian(8)
```

```julia
function solve(param)
    a, b, N, α, β, f = param
    x = range(a, b; length=N+1)
    h = step(x)
    A = Laplacian(N)
    b = BC(N, h, α, β)
    v = f.(x[2:end])
    u = [α; A\(h^2*v - b)] # solve h²(A + b)u = v
    (; x, u, param)
end
```

```julia
function plot_sol(sol, u_exact; kwargs...)
    x, u, param = sol
    a, b, N, α, β, f = param
    xs = range(a, b; length=1001)
    
    P = plot(x, u; label="numerical")
    isnothing(u_exact) || plot!(xs, u_exact.(Ref(param), xs); label="exact", ls=:dash)
    title!("N = $N"; titlefontsize=10)
    
    fs = f.(xs)
    ylim = extrema(fs)
    ydiff = ylim[2] - ylim[1]
    ylim =  (ylim[1] - 0.1(ydiff + 1), ylim[2] + 0.1(ydiff + 1))
    title_str = @sprintf("u(%.2f) = %.2f,  u'(%.2f) = %.2f", a, α, b, β) 
    Q = plot(xs, fs; label="f(x)", ylim)
    title!(title_str; titlefontsize=10)
    
    plot(P, Q; size=(800, 300))
end
```

```julia
param1a = (
    a = -1.0,
    b = 1.0,
    N = 3,
    α = 1.0,
    β = -1.0,
    f = x -> 1.0,
)

"""Assume f is a constant function.""" 
function u1_exact(param, x)
    a, b, N, α, β, f = param
    c = f(0.0)
    (c/2)x^2 + (-c*b + β)x - (c/2)a^2 - (-c*b + β)a + α
end

sol1a = solve(param1a)
plot_sol(sol1a, u1_exact)
```

```julia
param1b = (
    a = -1.0,
    b = 1.0,
    N = 3,
    α = 1.0,
    β = 0.9,
    f = x -> 1.0,
)

sol1b = solve(param1b)
plot_sol(sol1b, u1_exact)
```

```julia
param2 = (
    a = 0.0,
    b = 2π,
    N = 10,
    α = -1.0,
    β = -1.5,
    f = sin,
)

"""Assume a = 0,  b = 2π,  f = sin"""
function u2_exact(param, x)
    a, b, N, α, β, f = param
    α + (β + 1)x - sin(x)
end

sol2 = solve(param2)
plot_sol(sol2, u2_exact)
```

```julia

```
