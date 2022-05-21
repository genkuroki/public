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

* https://cjdoris.github.io/Bokeh.jl/stable/
* https://github.com/cjdoris/Bokeh.jl

```julia
using Bokeh
n = 2_000
z = rand(1:3, n)
x = randn(n) .+ [-2, 0, 2][z]
y = randn(n) .+ [-1, 3, -1][z]
color = Bokeh.PALETTES["Julia3"][z]
p = figure(title="Julia Logo")
plot!(p, Scatter; x, y, color, alpha=0.4, size=10)
p
```

```julia

```
