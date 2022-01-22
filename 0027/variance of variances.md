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
    display_name: Julia 1.7.1
    language: julia
    name: julia-1.7
---

```julia
using Random
using Distributions
using Plots
plot(sin)
```

```julia
function diffvarvar(dist, n)
    σ = std(dist)
    μ₄ = kurtosis(dist) + 3
    σ^4/n^2 * ((2 - 1/n)*μ₄ - (3 - (5n-3)/(n*(n-1))))
end

function diffvarvarnormal(dist, n)
    σ = std(dist)
    (3n - 1)/(n^2*(n-1)) * σ^4
end

function var_unbiasedvar(dist, n)
    σ = std(dist)
    μ₄ = kurtosis(dist) + 3
    σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
end

function var_var(dist, n)
    σ = std(dist)
    ((n-1)/n)^2 * var_unbiasedvar(dist, n) + σ^4/n^2
end

function simvar(dist, n; L = 10^7)
    σ² = var(dist)
    X = Vector{Float64}(undef, n)
    V = [var(rand!(dist, X); corrected=false) for _ in 1:L]
    U = [var(rand!(dist, X)) for _ in 1:L]
    varU = var(U; mean=σ²)
    varV = var(V; mean=σ²) 
    varU, varV, varU - varV
end
```

```julia
dist = Normal(1, 2)
n = 10
simvar(dist, n; L = 10^8), var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), diffvarvarnormal(dist, n)
```

```julia
dist = Exponential(2)
n = 10
simvar(dist, n; L = 10^8), var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), diffvarvarnormal(dist, n)
```

```julia
using SymPy
@vars σ μ₄ n

σ = 1
μ₄ = 0
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
([vu, vv, vu - vv] .|> simplify) .* -n^2 .|> simplify
```

```julia

```
