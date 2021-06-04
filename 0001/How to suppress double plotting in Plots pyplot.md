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
using Plots
pyplot()

if backend() == Plots.PyPlotBackend() && !@isdefined(IJulia_Setting_Changed)
    IJulia.pop_postexecute_hook(PyPlot.display_figs)
    IJulia.push_postexecute_hook(PyPlot.close_figs)
    IJulia_Setting_Changed = nothing
end
@show IJulia.postexecute_hooks

using MyUtils: showimg
```

```julia
plot(sin; label="", title="y = sin(x)", titlefontsize=16)
```

```julia
plot(sin; label="", title="y = sin(x)", titlefontsize=16, dpi=600)
savefig("Untitled.png")
showimg("image/png", "Untitled.png"; tag="img width=\"600\" height=\"400\"")
```

```julia
f(x, y) = exp(-x^2 + 1.5x*y - y^2)
x = y = range(-3, 3; length=300)
heatmap(x, y, f; size=(500, 400))
```

```julia
heatmap(x, y, f; size=(500, 400), dpi=600)
savefig("Untitled.png")
showimg("image/png", "Untitled.png"; tag="img width=\"500\" height=\"400\"")
```

```julia

```
