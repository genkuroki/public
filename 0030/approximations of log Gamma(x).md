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

```julia
using Plots
default(fmt = :png)
using SpecialFunctions
```

```julia
f(x) = loggamma(x)
f0(x) = x * log(x) - x - (1/2)*log(x) + log(√(2π))
f1(x) = x * log(x) - x - (1/2)*log(x)
f2(x) = x * log(x) - x
f3(x) = x * log(x)

function plot_loggamma(N)
    x = range(1, N, 1000)
    plot(; legend=:outertop)
    plot!(x, f; label="y = log Γ(x)")
    plot!(x, f0; label="y = x log(x) - x - (1/2)log(x) + log(√(2π))", ls=:dash)
    plot!(x, f1; label="y = x log(x) - x - (1/2)log(x)", ls=:dot, lw=1.5)
    plot!(x, f2; label="y = x log(x) - x", ls=:dashdot)
    plot!(x, f3; label="y = x log(x)")
end
```

```julia
plot_loggamma(3)
```

```julia
plot_loggamma(10)
```

```julia
plot_loggamma(30)
```

```julia
plot_loggamma(6.022e23)
```

```julia

```
