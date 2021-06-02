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
`Laplacian(N)` is the 1D discrete Laplacian
with left Dirichlet and right Neumann boundary conditions.

The 1D Poisson equation
\$u''(x) + f(x) = 0\$, 
\$u(a) = \\alpha\$,
\$u'(b) = \\beta\$
is approximated by \$h^2(A u + b) = v\$,
where
\$x = \\mathrm{range}(a, b; \\mathrm{length}=N+1)\$,
\$h = \\mathrm{step}(x)\$,
\$A = \\mathrm{Laplacian}(N)\$, 
\$b = \\mathrm{BC}(N, h, \\alpha, \\beta)\$, 
and
\$v = f.(x[2{:}\\mathrm{end}])\$.

__Remark.__ Laplacian(N) × (-1) is equal to the Cartan matrix of B-type (or C-type).
"""
Laplacian(N) = Tridiagonal(
    [fill(1, N-2); 2],
    [fill(-2, N-1); -2],
    fill(1, N-1)
)

BC(N, h, α, β) = [α; zeros(N-2); 2β*h]
```

```julia
?Laplacian
```

```julia
Laplacian(8)
```

```julia
"""
`solve_1d_poisson(param)` solves the discrete version of the 1D Poisson equation
\$u''(x) + f(x) = 0\$, \$u(a) = \\alpha\$, \$u'(b) = \\beta\$,
where `param = (a, b, α, β, f, N)`
"""
function solve_1d_poisson(param)
    a, b, α, β, f, N = param
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
?solve_1d_poisson
```

```julia
"""
`plot_1d_poisson(sol, u_exact)` plots the numerical solution `sol` and the exact solution function `u_exact`.
"""
function plot_1d_poisson(sol, u_exact)
    x, u, param = sol
    a, b, α, β, f, N = param
    xs = range(a, b; length=1001)
    
    P = plot(x, u; label="numerical")
    isnothing(u_exact) || plot!(xs, u_exact.(Ref(param), xs); label="exact", ls=:dash)
    title_str = @sprintf("u(%.2f) = %.2f,  u'(%.2f) = %.2f,  N = %d", N , a, α, b, β) 
    title!(title_str; titlefontsize=10)
    
    fs = f.(xs)
    ylim = extrema(fs)
    ydiff = ylim[2] - ylim[1]
    ylim =  (ylim[1] - 0.1(ydiff + 1), ylim[2] + 0.1(ydiff + 1))
    Q = plot(xs, fs; label="y = f(x)", ylim)
    title!("u''(x) = f(x)"; titlefontsize=10)
    
    plot(P, Q; size=(800, 300))
end
```

```julia
?plot_1d_poisson
```

```julia
param1a = (
    a = -1.0,
    b = 1.0,
    α = 1.0,
    β = -1.0,
    f = x -> 1.0,
    N = 3,
)

"""Assume f is a constant function.""" 
function u1_exact(param, x)
    a, b, α, β, f, N = param
    c = f(0.0)
    (c/2)x^2 + (-c*b + β)x - (c/2)a^2 - (-c*b + β)a + α
end

sol1a = solve_1d_poisson(param1a)
plot_1d_poisson(sol1a, u1_exact)
```

```julia
param1b = (
    a = -1.0,
    b = 1.0,
    α = 1.0,
    β = 0.9,
    f = x -> 1.0,
    N = 3,
)

sol1b = solve_1d_poisson(param1b)
plot_1d_poisson(sol1b, u1_exact)
```

```julia
param2 = (
    a = 0.0,
    b = 2π,
    α = -1.0,
    β = -1.5,
    f = sin,
    N = 10,
)

"""Assume a = 0,  b = 2π,  f = sin"""
function u2_exact(param, x)
    a, b, α, β, f, N = param
    α + (β + 1)x - sin(x)
end

sol2 = solve_1d_poisson(param2)
plot_1d_poisson(sol2, u2_exact)
```

```julia

```
