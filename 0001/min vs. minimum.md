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

https://twitter.com/julia_kizi/status/1400703388363468801

```julia
using BenchmarkHistograms
```

```julia
A = randn(10^5)

f(A) = min(A...)
g(A) = minimum(A)
```

```julia
@benchmark f($A)
```

```julia
@benchmark g($A)
```

```julia

```
