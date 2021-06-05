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
?gcdx
```

```julia
xgcd(a) = (abs(a), sign(a))
xgcd(a, b) = gcdx(a, b)
function xgcd(a, b, c...)
    d, u, v = xgcd(a, b)
    X = xgcd(d, c...)
    g, w, x = X[1], X[2], X[3:end]
    (g, w*u, w*v, x...)
end
```

```julia
A = [2310*rand(1:10^5), 2310*rand(1:10^5), 210*rand(1:10^5), 30*rand(1:10^5), 6*rand(1:10^5)]
@show A xgcd(A...); 
```

```julia
collect(xgcd(A...)[2:end])'A
```

```julia

```
