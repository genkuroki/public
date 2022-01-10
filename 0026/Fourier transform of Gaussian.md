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
using SymPy: SymPy, sympy, @syms, PI, oo
@syms x::real y::real ξ::real
```

```julia
sympy.integrate(exp(-2π*im*x*ξ - π*x^2), (x, -oo, oo)).simplify()
```

```julia
sympy.integrate(exp(-2PI*im*x*ξ - π*x^2), (x, -oo, oo)).simplify()
```

```julia
sympy.integrate(exp(y^2), (y, x, 1)).simplify()
```

```julia
sympy.integrate(exp(y^2), (y, x, 1), (x, 0, 1)).simplify()
```

https://github.com/genkuroki/public/blob/d54bc52c4198e88468e85233130882222890c02f/0023/Maxima.jl%20example.ipynb

```julia
using Maxima: Maxima, MExpr, @m_str, mcall
Base.convert(::Type{Any}, x::MExpr) = x
```

```julia
m"integrate(integrate(exp(y^2), y, x, 1), x, 0, 1)" |> mcall
```

```julia
m"integrate(integrate(exp(y^2), x, 0, y), y, 0, 1)" |> mcall
```

```julia
m"integrate(log(sin(x)), x)" |> mcall
```

```julia

```
