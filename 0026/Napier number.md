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
    display_name: Julia 1.7.0
    language: julia
    name: julia-1.7
---

```julia
f(n) = (1 + 1/n)^n
[(n, f(n), exp(1)) for n in 10 .^ (1:7)]
```

```julia
[(n, f(n)/exp(1) - 1, 1/(2n)) for n in 10 .^ (1:7)]
```

```julia

```
