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
(1:9)' .* (1:9)
```

```julia
using Plots
f(x, y) = exp(-x^2 + x*y/2 - y^2)
x = y = range(-3, 3; length=301)
z = f.(x', y)
surface(x, y, z; colorbar=false, size=(720, 540), camera=(60, 60), color=:CMRmap)
```

```julia
using Plots
f(x, y) = exp(-x^2 + x*y/2 - y^2)
x = y = range(-3, 3; length=301)
z = f.(x', y)
surface(x, y, z; colorbar=false, size=(720, 540), camera=(60, 60), color=:gist_earth)
```

```julia

```
