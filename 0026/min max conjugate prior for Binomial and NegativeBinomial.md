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
using Memoization

default(size=(500, 400))
plot(sin);
```

```julia
@memoize function E_bin(f, n, p)
    bin = Binomial(n, p)
    sum(f(k) * pdf(bin, k) for k in support(bin))
end

@memoize function E_negbin(g, k, p)
    negbin = LocationScale(k, 1, NegativeBinomial(k, p))
    m, s = mean(negbin), std(negbin)
    nmax = round(Int, m + 5s)
    sum(g(n) * pdf(negbin, n) for n in k:nmax)
end

function kl(n, k, p; a = 0.5, b = 0.5)
    p̂ = (k + a)/(n + a + b)
    -(p*log(p̂) + (1 - p)*log(1 - p̂)) - entropy(Bernoulli(p)) 
end

function squared_error(n, k, p; a = 0.5, b = 0.5)
    p̂ = (k + a)/(n + a + b)
    (p̂ - p)^2
end

function abs_error(n, k, p; a = 0.5, b = 0.5)
    p̂ = (k + a)/(n + a + b)
    abs(p̂ - p)
end
```

```julia
n = 10
@memoize f(n, a, b) =  maximum(p -> E_bin(k -> kl(n, k, p; a, b), n, p), 0.0:0.01:1)
a = b = 0.01:0.01:1
@time z = f.(n, a, b')
@show val, idx = findmin(z)
@show a[idx[1]], b[idx[2]]
P = contourf(a, b, z'; label="")
Q = surface(a, b, z'; label="")
plot(P, Q; size=(800, 400), colorbar=false, camera=(60, 60))
```

```julia
k = 3
@memoize g(k, a, b) =  maximum(p -> E_negbin(n -> kl(n, k, p; a, b), k, p), 0.01:0.01:0.99)
a = 0.2:0.002:0.5
b = 0.3:0.002:0.55
@time z = g.(k, a, b')
@show val, idx = findmin(z)
@show a[idx[1]], b[idx[2]]
P = contourf(a, b, z'; label="")
Q = surface(a, b, z'; label="")
plot(P, Q; size=(800, 400), colorbar=false, camera=(60, 60))
```

```julia
n = 10
@memoize F(n, a, b) =  maximum(p -> E_bin(k -> squared_error(n, k, p; a, b), n, p), 0:0.01:1)
a = b = 0.5:0.01:2
@time z = F.(n, a, b')
@show val, idx = findmin(z)
@show a[idx[1]], b[idx[2]]
P = contourf(a, b, z'; label="")
Q = surface(a, b, z'; label="")
plot(P, Q; size=(800, 400), colorbar=false, camera=(60, 60))
```

```julia
k = 3
@memoize G(k, a, b) =  maximum(p -> E_negbin(n -> squared_error(n, k, p; a, b), k, p), 0.01:0.01:0.99)
a = 0.01:0.002:0.25
b = 0.5:0.002:0.8
@time z = G.(k, a, b')
@show val, idx = findmin(z)
@show a[idx[1]], b[idx[2]]
P = contourf(a, b, z'; label="")
Q = surface(a, b, z'; label="")
plot(P, Q; size=(800, 400), colorbar=false, camera=(60, 60))
```

```julia
n = 10
@memoize f1(n, a, b) =  maximum(p -> E_bin(k -> abs_error(n, k, p; a, b), n, p), 0:0.01:1)
a = b = 0.5:0.01:2
@time z = f1.(n, a, b')
@show val, idx = findmin(z)
@show a[idx[1]], b[idx[2]]
P = contourf(a, b, z'; label="")
Q = surface(a, b, z'; label="")
plot(P, Q; size=(800, 400), colorbar=false, camera=(60, 60))
```

```julia
k = 3
@memoize g1(k, a, b) =  maximum(p -> E_negbin(n -> abs_error(n, k, p; a, b), k, p), 0.01:0.01:0.99)
a = 0.2:0.01:1
b = 0.5:0.001:0.7
@time z = g1.(k, a, b')
@show val, idx = findmin(z)
@show a[idx[1]], b[idx[2]]
P = contourf(a, b, z'; label="")
Q = surface(a, b, z'; label="")
plot(P, Q; size=(800, 400), colorbar=false, camera=(60, 60))
```

```julia

```
