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

https://medium.com/@andreaskuhn92/how-to-solve-the-same-numerical-problem-in-7-different-programming-languages-a64daac3ed64

```julia
using BenchmarkTools
```

```julia
function julia_logistic(N,dt,u)
    @fastmath @inbounds begin 
        # Parameters
        u0 = 1e-5
        # Right hand side function
        f(U)= U*(1.0-U)
        # Discretization
        u[1] = u0
        for k = 1:(N-1)
            u[k+1] = u[k] + dt*f(u[k])
        end
        return(u)
    end
end

N = 10^3
dt = 25/N
u = zeros(N)
u = @btime julia_logistic(N, dt, $u)
```

```julia
using LoopVectorization

function julia_logistic_turbo(N,dt,u)
    # Parameters
    u0 = 1e-5
    # Right hand side function
    f(U)= U*(1.0-U)
    # Discretization
    u[1] = u0
    @turbo for k = 1:(N-1)
        u[k+1] = u[k] + dt*f(u[k])
    end
    return(u)
end

N = 10^3
dt = 25/N
v = zeros(N)
v = @btime julia_logistic_turbo(N, dt, $v)
```

```julia
u ≈ v
```

```julia

```
