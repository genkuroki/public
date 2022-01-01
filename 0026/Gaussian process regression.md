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
using LinearAlgebra
using Distributions
using Plots
default(fmt=:png, ms=3, titlefontsize=10, tickfontsize=6)
plot(sin; size=(200, 150));
```

## 単なるGauss過程 = データ無しの段階でのGauss過程回帰

```julia
mu(x) = zero(x)
#mu(x) = sin(x)
g(x, x′; a=1, b=1, σ=0.2) = a * exp(-(x - x′)^2/(2b^2)) + σ^2 * (x == x′)

L = 600
x = range(-0.1π, 6.1π, L)
μ = mu.(x)
Σ = Symmetric(@.(g(x, x')))
mvn = MvNormal(μ, Σ)
```

```julia
y = rand(mvn)
plot(x, y; label="")
```

```julia
nsamples = 10
y = rand(mvn, nsamples)
plot(; title="$nsamples samples")
plot!(x, y; label="", lw=0.5, alpha=0.5)
```

## データから得られる条件付き確率分布の構成 = Gauss過程回帰

```julia
# テストデータの生成

n = 10
X = range(0, 6π, n)
Y = sin.(X) + 0.3randn(n)

L = 600
x = range(-0.1π, 6.1π, L) # 予測したい x の範囲

plot(; title="test data")
scatter!(X, Y; label = "", c=:red)
plot!(x, sin.(x); label="", c=:orange)
```

```julia
# Gauss過程回帰 = データ (X, Y) から得られる条件付き確率分布の構成

μ_X, μ_x = mu.(X), mu.(x)
Σ_XX, Σ_Xx, Σ_xX, Σ_xx = g.(X, X'), g.(X, x'), g.(x, X'), g.(x, x')
A = Σ_xX / Σ_XX
μ = μ_x + A * (Y - μ_X)
Σ = Symmetric(Σ_xx - A * Σ_Xx)
mvn_pred = MvNormal(μ, Σ)
```

```julia
y = rand(mvn_pred)

plot(; legend=:outertop)
plot!(x, y; label="", title="sample of predictive distribution")
scatter!(X, Y; label = "", c=:red)
plot!(sin, extrema(X)...; label="", c=:orange)
```

```julia
nsamples = 20
y = rand(mvn_pred, nsamples)

plot(; title="$nsamples samples of predictive distribution")
plot!(x, y; label="", lw=0.3, alpha=0.5)
scatter!(X, Y; label = "", c=:red)
plot!(sin, extrema(X)...; label="", c=:orange)
```

```julia
s = diag(Σ)
plot(; title="μ ± 2σ interval of predictive distribution")
plot!(x, [μ μ]; label="", fillrange=[μ-2s μ+2s], c=1, fillalpha=0.2)
scatter!(X, Y; label = "", c=:red)
plot!(sin, extrema(X)...; label="", c=:orange)
```

```julia

```
