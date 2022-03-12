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
using Distributions
using StatsPlots
default(fmt=:png)
using Random
```

```julia tags=[]
function biased_var_vs_unbiased_var(; dist = Normal(), n = 10, L = 1000,
        loss_func = (x, v) -> (x - v)^2)
    σ² = var(dist)
    Loss_biased = Vector{Float64}(undef, L)
    Loss_unbiased = Vector{Float64}(undef, L)
    tmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        sample = rand!(dist, tmp[Threads.threadid()])
        b = var(sample; corrected=false) # biased var
        u = var(sample; corrected=true)  # unbiased var
        Loss_biased[i] = loss_func(b, σ²)
        Loss_unbiased[i] = loss_func(u, σ²)
    end
    Loss_biased, Loss_unbiased
end

function plot_score(; dist = Normal(), n = 10, L = 10^4,
        loss_func = (x, v) -> (x - v)^2,
        bin = 100)
    Loss_biased, Loss_unbiased = biased_var_vs_unbiased_var(; dist, n, L, loss_func)
    Score = Loss_unbiased - Loss_biased
    @show mean(Score) std(Score)
    P1 = plot(cumsum(Score); label="total score", legend=:outertop)
    P2 = histogram(Score; alpha=0.3, label="score", legend=:outertop, bin)
    vline!([0]; label="", c=:red)
    plot(P1, P2; size=(800, 300))
end

Random.seed!(4649373)
```

```julia
plot_score(; loss_func = (x, v) -> abs(x - v), L = 100, bin = 20)
```

```julia
plot_score(; loss_func = (x, v) -> abs(x - v))
```

```julia
plot_score(; loss_func = (x, v) -> (x - v)^2, L = 100, bin = 20)
```

```julia
plot_score(; loss_func = (x, v) -> (x - v)^2)
```

```julia
@time plot_score(; loss_func = (x, v) -> (x - v)^2, L = 10^8)
```

```julia
plot_score(; loss_func = (x, v) -> abs(√x - √v))
```

```julia
plot_score(; loss_func = (x, v) -> (√x - √v)^2)
```

```julia
plot_score(; loss_func = (x, v) -> abs(log(x) - log(v)))
```

```julia
plot_score(; loss_func = (x, v) -> (log(x) - log(v))^2)
```

```julia
plot_score(; loss_func = (x, v) -> abs(log(x) - log(v)), n = 100)
```

```julia
plot_score(; loss_func = (x, v) -> (log(x) - log(v))^2, n = 100)
```

```julia

```
