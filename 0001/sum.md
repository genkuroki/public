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

```julia
f(x, y) = x * y
X = Y = 1:9
```

```julia
sum(f(x, y) for x in X for y in Y)
```

```julia
sum(f(x, y) for x in X, y in Y)
```

```julia
sum(t -> f(t...), Iterators.product(X, Y))
```

```julia
sum(x -> x^3, X)
```

```julia

```
