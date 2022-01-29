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

# 不偏分散の分散と標本分散の平均二乗誤差

* 黒木玄
* 2022-01-24

$
\newcommand\var{\operatorname{var}}
$

$X_1, X_2, \ldots, X_n$ はその各々が分散 $\sigma^2$ と尖度 $\kappa$ (正規分布で0になるようにしたもの, kurtosis)を持つ分布に従う独立同分布確率変数達であるとする。このとき, その標本平均 $\bar{X}$, 不偏分散 $U$ と標本分散 $V$ が次のように定義される:

$$
\bar{X} = \frac{1}{n}\sum_{i=1}^n X_i, \quad
U = \frac{1}{n-1}\sum_{i=1}^n (X_i - \bar{X})^2, \quad
V = \frac{1}{n}\sum_{i=1}^n (X_i - \bar{X})^2 = \frac{n-1}{n}U.
$$


このとき, $\var(\bar{X}) = E[(\bar{X}-\mu)^2] = {\sigma^2}/{n}$ でかつ, 

$$
E[U] = \sigma^2, \quad
\var(U) = \sigma^4\left(\frac{\kappa+3}{n} - \frac{n-3}{n(n-1)}\right), \quad
E\left[\left(V - \sigma^2\right)^2\right] = \left(\frac{n-1}{n}\right)^2 \var(U) + \frac{\sigma^4}{n^2}.
$$

ゆえに

$$
\var(U) - E\left[\left(V - \sigma^2\right)^2\right] =
\frac{\sigma^4}{n^2}
\left(\left(2 - \frac{1}{n}\right)(\kappa + 3) - \left(3 - \frac{5n-3}{n(n-1)}\right)\right).
$$


これより, $\kappa > -\dfrac{3}{2}$ ならば十分大きな $n$ について $\var(U) > E\left[\left(V - \sigma^2\right)^2\right]$.


正規分布の場合には $\kappa = 0$ となるので,

$$
\var(U) = \frac{2\sigma^4}{n-1}, \quad
\var(U) - E\left[\left(V - \sigma^2\right)^2\right] = \frac{(3n - 1)\sigma^4}{n^2(n-1)} > 0.
$$

```julia
using Random
using Distributions
using StatsPlots
plot(TDist(1), -6, 6; size=(200, 150))
```

```julia
function diffvarvar(dist, n)
    σ = std(dist)
    μ₄ = kurtosis(dist) + 3
    σ^4/n^2 * ((2 - 1/n)*μ₄ - (3 - (5n-3)/(n*(n-1))))
end

function diffvarvarnormal(dist, n)
    σ = std(dist)
    (3n - 1)/(n^2*(n-1)) * σ^4
end

function var_unbiasedvar(dist, n)
    σ = std(dist)
    μ₄ = kurtosis(dist) + 3
    σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
end

function var_var(dist, n) # 正確には分散ではなく、σ²に対する平均二乗誤差
    σ = std(dist)
    ((n-1)/n)^2 * var_unbiasedvar(dist, n) + σ^4/n^2
end

function simvar(dist, n; L = 10^7)
    σ² = var(dist)
    tmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    V = Vector{Float64}(undef, L)
    U = similar(V)
    Threads.@threads for i in 1:L
        X = rand!(dist, tmp[Threads.threadid()])
        V[i] = var(X; corrected=false)
        U[i] = var(X)
    end
    varU = var(U; mean=σ²)
    varV = var(V; mean=σ²) 
    varU, varV, varU - varV, varV/varU, varU/varV, √(varU/varV)
end
```

```julia
dist = Normal()
n = 10
diffvarvarnormal(dist, n) |> display
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Normal(1, 2)
n = 10
diffvarvarnormal(dist, n) |> display
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Exponential()
n = 10
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Exponential(2)
n = 10
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Bernoulli(0.5)
n = 10
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Bernoulli(0.5)
n = 100
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Bernoulli(0.5)
n = 1000
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Bernoulli(0.4)
n = 10
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Bernoulli(0.3)
n = 10
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
dist = Bernoulli(0.2)
n = 10
simvar(dist, n; L = 10^8) |> display
var_unbiasedvar(dist, n), var_var(dist, n), diffvarvar(dist, n), var_var(dist, n)/var_unbiasedvar(dist, n), var_unbiasedvar(dist, n)/var_var(dist, n), √(var_unbiasedvar(dist, n)/var_var(dist, n))
```

```julia
using SymPy
@vars σ μ₄ n

σ = 1
μ₄ = 0
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
[vu, vv] .|> simplify .|> factor
```

```julia
(vu - vv)*(-n^2) |> simplify
```

```julia
@vars σ μ₄ n

μ₄ = 3 # normal dist. case
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
[vu, vv, vv/vu] .|> simplify .|> factor
```

```julia
@vars σ μ₄ n

μ₄ = 3 # normal dist. case
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
F(n => 10)
```

```julia
@vars σ μ₄ n

μ₄ = 1 # Bernoulli(1/2)
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
F(n => 10)
```

```julia
@vars σ μ₄ n

μ₄ = 3//2
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
@vars σ μ₄ n

μ₄ = 0
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
@vars σ μ₄ n

μ₄ = 1//2
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
@vars σ μ₄ n

μ₄ = 1
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
@vars σ μ₄ n

μ₄ = 3//2
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
@vars σ μ₄ n

μ₄ = 2
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
@vars σ μ₄ n

μ₄ = 5//2
vu = σ^4 * (μ₄/n - (n-3)/(n*(n-1)))
vv =  ((n-1)/n)^2 * vu + σ^4/n^2
F = vv/vu |> simplify |> factor
```

```julia
skewness(Normal()), skewness(Gamma(20, 1)), skewness(InverseGamma(20, 1))
```

```julia
kurtosis(Normal()), kurtosis(Gamma(20, 1)), kurtosis(InverseGamma(20, 1)), kurtosis(Laplace())
```

```julia
kurtosis(TDist(4.0001)), kurtosis(TDist(10)), kurtosis(TDist(20))
```

```julia
skewness(Beta(0.1, 9.9)), skewness(Beta(7, 3)), skewness(Beta(7, 3)), skewness(Beta(9.9, 0.1))
```

```julia
kurtosis(Beta(0.1, 9.9)), kurtosis(Beta(7, 3)), kurtosis(Beta(7, 3)), kurtosis(Beta(9.9, 0.1))
```

```julia
kurtosis.(Bernoulli.(0.1:0.1:0.9))
```

```julia
Binomial.(100, [0.01, 0.3, 0.7, 0.99]) .|> skewness
```

```julia
Binomial.(100, [0.01, 0.3, 0.7, 0.99]) .|> kurtosis
```

```julia

```
