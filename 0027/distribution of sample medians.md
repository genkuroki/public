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
using Distributions
using StatsPlots
default(fmt=:png)
using Random

mediannormal(dist, n) = Normal(median(dist), 1/√(4n*pdf(dist, median(dist))^2))

function plot_mediandist(; dist = Uniform(), n = 100, L = 10^5)
    Median = Vector{Float64}(undef, L)
    tmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        Median[i] = median(rand!(dist, tmp[Threads.threadid()]))
    end
    @show dist
    @show n
    @show mean(Median) median(dist)
    @show var(Median) 1/(4n*pdf(dist, median(dist))^2)
    histogram(Median; norm=true, alpha=0.3, label="dist. of sample medians")
    plot!(mediannormal(dist, n); ls=:dash, lw=1.5, label="normal approx.")
end
```

```julia
plot_mediandist(; dist = Uniform(), n = 100, L = 10^5)
```

```julia
plot_mediandist(; dist = Beta(0.5, 0.6), n = 100, L = 10^5)
```

```julia
plot_mediandist(; dist = Gamma(2, 1), n = 100, L = 10^5)
```

```julia
plot_mediandist(; dist = LogNormal(), n = 100, L = 10^5)
```

```julia

```
