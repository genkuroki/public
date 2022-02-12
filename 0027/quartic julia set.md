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

* https://twitter.com/1Hassium/status/1490623409637949442
* https://twitter.com/1Hassium/status/1492019379021975555

```julia
using Plots
default(legend=false, colorbar=false, titlefontsize=12)
gr(fmt=:png)
using CUDA
```

```julia
plotjulia(j; size=(540, 540), color=:gist_earth, kwargs...) = heatmap(j; color, size,
    colorbar=false, ticks=false, frame=false, axis=false, bgcolor=:black, kwargs...)

plotmandelbrot(m; size=(540, 540), color=reverse(cgrad(:jet1)), kwargs...) = plotjulia(m; color, size, kwargs...)

function julia(f, z, c, maxiters=2^10, threshold2=Inf, atol2=eps())
    w = z
    for i in 1:maxiters+1
        w = f(w, c)
    end
    for i in 1:maxiters
        z = f(z, c)
        abs2(z) ≥ threshold2 && return float(i)
        abs2(z - w) < atol2 && return float(i)
    end
    NaN
end

mandelbrot(f, maxiters=2^10, threshold2=Inf, atol2=eps()) =
    julia(f, zero(c), c, maxiters, threshold2, atol2)
```

```julia
n = 2^9
x = range(-1.5, 1.5; length=n)
y = range(-1.5, 1.5; length=n)
z = complex.(x', y)
c = -0.380 + 0.605im
f_julia(z, c) = z^2 + c

@show typeof(z)
@time j = julia.(f_julia, z, c, 2^8)
@time j = julia.(f_julia, z, c, 2^8)
@show typeof(j)
plotjulia(j)
```

```julia
g(z, c) = z^4 + c
c = 0.63 + 0.37im
```

```julia
x = range(-1, 1, 201)
y = range(-1, 1, 201)
z = complex.(x', y)
@show typeof(z)
@time j = julia.(g, z, c, 2^14, 1e2, 1e-5)
@time j = julia.(g, z, c, 2^14, 1e2, 1e-5)
@show typeof(j)
m = j
plotmandelbrot(.∛m; color=:jet1)
```

```julia
x = range(-1, 1, 201)
y = range(-1, 1, 201)
z = complex.(x', y)
@show typeof(z)
@time z_cuda = CuMatrix{ComplexF64}(z)
@show typeof(z_cuda)
@time j_cuda = julia.(g, z_cuda, c, 2^14, 1e2, 1e-5)
@time j_cuda = julia.(g, z_cuda, c, 2^14, 1e2, 1e-5)
@show typeof(j_cuda)
@time m = Matrix(j_cuda)
plotmandelbrot(.∛m; color=:jet1)
```

```julia
x = range(-1, 1, 2001)
y = range(-1, 1, 2001)
z = complex.(x', y)
@show typeof(z)
@time z_cuda = CuMatrix{ComplexF64}(z)
@show typeof(z_cuda)
@time j_cuda = julia.(g, z_cuda, c, 2^14, 1e2, 1e-5)
@show typeof(j_cuda)
@time m = Matrix(j_cuda)
plotmandelbrot(.∛m; color=:jet1)
```

```julia
plotmandelbrot(mod.(m, 2^11); color=:jet1)
```

```julia

```
