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
using Distributions
using Optim
using Plots
default(size=(400, 300), titlefontsize=10, tickfontsize=6)
plot(sin; size=(200, 150))
```

```julia
function f(X, m, w)
    μ = [w[1], w[2]]
    Σ = [w[3] w[4]; w[4] w[5]]
    Y = copy(X)
    Y[2, 1:m] = w[6:5+m]
    -loglikelihood(MvNormal(μ, Σ), Y)
end

dist = MvNormal([0, 0], [2 -1; -1 2])
n = 2^7
X = rand(dist, n)

m = 2^6
o = optimize(w -> f(X, m, w), [0;0; 1;0;1; zeros(m)], Optim.Options(iterations=10^6))
```

```julia
w = o.minimizer

Y = copy(X)
Y[2, 1:m] = w[6:5+m]

@show w[1:2] w[3:5];
```

```julia
fit_mle(MvNormal, X)
```

```julia
fit_mle(MvNormal, X[:, m+1:end])
```

```julia
fit_mle(MvNormal, Y)
```

```julia
P = scatter(X[1,:], X[2,:])

Q = scatter(Y[1,21:end], Y[2,21:end])
scatter!(Y[1,1:20], Y[2,1:20])

plot(P, Q; size=(800, 300))
```

```julia

```
