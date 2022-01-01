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
using StatsPlots
default(size=(400, 300), titlefontsize=10, tickfontsize=6)
plot(sin);
```

```julia
f(n, p, k) = cdf(Binomial(n, p), k)
g(n, p, k) = ccdf(Beta(k+1, n-k), p)

n = 20
k = 7
p = range(0, 1, 1000)
plot(p, f.(n, p, k))
plot!(p, g.(n, p, k), ls=:dash)
```

```julia
F(n, p, k) = ccdf(NegativeBinomial(k, p), n-k-1)
G(n, p, k) = ccdf(Beta(k, n-k), p)

n = 20
k = 7
p = range(eps(), 1, 1000)
plot(p, F.(n, p, k))
plot!(p, G.(n, p, k), ls=:dash)
```

```julia
λ = 7
N = 100
poi = Poisson(λ)
bin = Binomial(N, λ/N)
plot(x -> cdf(poi, x), 0, 15; label="")
plot!(x -> cdf(bin, x), 0, 15; label="", ls=:dash)
```

```julia
λ = 7.89
expon = Exponential(1/λ)

N = 10^2
geom = LocationScale(1/N, 1/N, Geometric(λ/N))

plot(x -> cdf(expon, x), 0, 1; label="")
plot!(x -> cdf(geom, x), 0, 1; label="", ls=:dash)
```

```julia
α = 12.345
λ = 7.89
gamma = Gamma(α, 1/λ)

N = 10^2
negbin = LocationScale(α/N, 1/N, NegativeBinomial(α, λ/N))

plot(x -> cdf(gamma, x), 0, 4; label="")
plot!(x -> cdf(negbin, x), 0, 4; label="", ls=:dash)
```

```julia
f(λ, t, n) = cdf(Poisson(λ*t), n-1)
g(λ, t, n) = ccdf(Gamma(n, 1/λ), t)

λ = 7
n = 10
plot(t -> f(λ, t, n), 0, 4; label="")
plot!(t -> g(λ, t, n), 0, 4; label="", ls=:dash)
```

```julia
f(λ, t, n) = cdf(Poisson(λ*t), n-1)
g(λ, t, n) = ccdf(Gamma(n, 1/λ), t)

λ = 7
t = 0.5
plot(n -> f(λ, t, n), 1, 10; label="")
plot!(n -> g(λ, t, n), 1, 10; label="", ls=:dash)
```

```julia
# n回での成功回数がk以下
f(n, p, k) = cdf(Binomial(n, p), k)

# k+1回成功するまでn+1回以上
g(n, p, k) = ccdf(NegativeBinomial(k+1, p), (n+1)-(k+1)-1)

n = 10
k = 5
plot(p -> f(n, p, k), eps(), 1; label="")
plot!(p -> g(n, p, k), eps(), 1; label="", ls=:dash)
```

```julia
# N回での成功回数がn-1以下
F(N, p, n) = cdf(Binomial(N, p), n-1)

# n回成功するまでN+1回以上
G(N, p, n) = ccdf(NegativeBinomial(n, p), N+1-n-1)

N = 10
n = 5
plot(p -> F(N, p, n), eps(), 1; label="")
plot!(p -> G(N, p, n), eps(), 1; label="", ls=:dash)
```

```julia
# N回での成功回数がn-1以下
F(N, p, n) = cdf(Binomial(N, p), n-1)

# n回成功するまでN+1回以上
G(N, p, n) = ccdf(NegativeBinomial(n, p), N+1-n-1)

N = 100
p = 0.1
plot(n -> F(N, p, n), 1, 20; label="")
plot!(n -> G(N, p, n), 1, 20; label="", ls=:dash)
```

```julia
f(λ, t, n) = cdf(Poisson(λ*t), n-1)
g(λ, t, n) = ccdf(Gamma(n, 1/λ), t)

# N回での成功回数がn-1以下
F(N, p, n) = cdf(Binomial(N, p), n-1)

# n回成功するまでN+1回以上
G(N, p, n) = ccdf(NegativeBinomial(n, p), N+1-n-1)

t = 1.23
λ = 6.78
N = 100
p = λ*t/N
plot(n -> F(N, p, n), 1, 20; label="")
plot!(n -> G(N, p, n), 1, 20; label="", ls=:dash)
plot!(n -> f(λ, t, n), 1, 20; label="", ls=:dashdot)
plot!(n -> g(λ, t, n), 1, 20; label="", ls=:dot)
```

```julia
f(λ, t, n) = cdf(Poisson(λ*t), n-1)
g(λ, t, n) = ccdf(Gamma(n, 1/λ), t)

# N回での成功回数がn-1以下
F(N, p, n) = cdf(Binomial(N, p), n-1)

# n回成功するまでN+1回以上
G(N, p, n) = ccdf(NegativeBinomial(n, p), N+1-n-1)

λ = 6.78
N = 30
n = 10
plot(t -> F(N, λ*t/N, n), 0, 3; label="")
plot!(t -> G(N, λ*t/N, n), eps(), 3; label="", ls=:dash)
plot!(t -> f(λ, t, n), 0, 3; label="", ls=:dashdot)
plot!(t -> g(λ, t, n), 0, 3; label="", ls=:dot)
```

```julia
f(λ, t, n) = cdf(Poisson(λ*t), n-1)
g(λ, t, n) = ccdf(Gamma(n, 1/λ), t)

# N回での成功回数がn-1以下
F(N, p, n) = cdf(Binomial(N, p), n-1)

# n回成功するまでN+1回以上
G(N, p, n) = ccdf(NegativeBinomial(n, p), N+1-n-1)

λ = 6.78
N = 10^2
n = 10
plot(t -> F(N, λ*t/N, n), 0, 3; label="")
plot!(t -> G(N, λ*t/N, n), eps(), 3; label="", ls=:dash)
plot!(t -> f(λ, t, n), 0, 3; label="", ls=:dashdot)
plot!(t -> g(λ, t, n), 0, 3; label="", ls=:dot)
```

```julia
f(λ, t, n) = cdf(Poisson(λ*t), n-1)
g(λ, t, n) = ccdf(Gamma(n, 1/λ), t)

# N回での成功回数がn-1以下
F(N, p, n) = cdf(Binomial(N, p), n-1)

# n回成功するまでN+1回以上
G(N, p, n) = ccdf(NegativeBinomial(n, p), N+1-n-1)

λ = 6.78
N = 10^3
n = 10
plot(t -> F(N, λ*t/N, n), 0, 3; label="")
plot!(t -> G(N, λ*t/N, n), eps(), 3; label="", ls=:dash)
plot!(t -> f(λ, t, n), 0, 3; label="", ls=:dashdot)
plot!(t -> g(λ, t, n), 0, 3; label="", ls=:dot)
```

```julia

```
