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
    display_name: Julia 1.8.0-beta1
    language: julia
    name: julia-1.8
---

```julia
using Distributions
using StatsPlots
default(fmt=:png)

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 20
L = 10^5
X = [mean(rand(dist, n)) for _ in 1:L]
histogram(X; norm=true, alpha=0.3, label="n=$n sample means")
plot!(Normal(μ, σ/√n); label="normal approx")
```

```julia
using Distributions
using StatsPlots

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 100
L = 10^5
X = [mean(rand(dist, n)) for _ in 1:L]
histogram(X; norm=true, alpha=0.3, label="n=$n sample means")
plot!(Normal(μ, σ/√n); label="normal approx")
```

```julia
plot(x -> pdf(dist, x), -3, 23; label="dist")
```

```julia
using Distributions
using StatsPlots
using StatsBase: ecdf
myecdf(A, x) = count(≤(x), A)/length(A)

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 20
L = 10^5
X = [mean(rand(dist, n)) for _ in 1:L]

f = ecdf(X)
g(x) = myecdf(X, x)
plot(; legend=:topleft)
plot!(x -> f(x), -2.5, 6.5; label="ecdf")
plot!(g; label="myecdf", ls=:dash)
plot!(x -> cdf(Normal(μ, σ/√n), x); label="normal approx", ls=:dashdot)
```

```julia
using Distributions
using StatsPlots
using Random

dist = Normal()
μ, σ = mean(dist), std(dist)
n = 20
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=2, msw=0, alpha=0.5, label="", title="Normal(), n=$n")
scatter!([(μ, σ^2)]; ms=4, label="")
```

```julia
using Distributions
using StatsPlots
using Random

dist = Uniform()
μ, σ = mean(dist), std(dist)
n = 20
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=2, msw=0, alpha=0.5, label="", title="Uniform(), n=$n")
scatter!([(μ, σ^2)]; ms=4, label="")
```

```julia
using Distributions
using StatsPlots
using Random

dist = Exponential()
μ, σ = mean(dist), std(dist)
n = 20
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=2, msw=0, alpha=0.5, label="", title="Exponential(), n=$n")
scatter!([(μ, σ^2)]; ms=4, label="")
```

```julia
using Distributions
using StatsPlots
using Random

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 20
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=1, msw=0, alpha=0.5, label="",
    title="MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05]), n=$n",
    titlefontsize=12)
scatter!([(μ, σ^2)]; ms=3, label="")
```

```julia
using Distributions
using StatsPlots
using Random

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 40
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=1, msw=0, alpha=0.5, label="",
    title="MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05]), n=$n",
    titlefontsize=12)
scatter!([(μ, σ^2)]; ms=3, label="")
plot!(; xlim=(-0.5, 4), ylim=(-2, 80))
```

```julia
using Distributions
using StatsPlots
using Random

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 80
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=1, msw=0, alpha=0.5, label="",
    title="MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05]), n=$n",
    titlefontsize=12)
scatter!([(μ, σ^2)]; ms=3, label="")
plot!(; xlim=(-0.5, 4), ylim=(-2, 80))
```

```julia
using Distributions
using StatsPlots
using Random

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 160
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=1, msw=0, alpha=0.5, label="",
    title="MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05]), n=$n",
    titlefontsize=12)
scatter!([(μ, σ^2)]; ms=3, label="")
plot!(; xlim=(-0.5, 4), ylim=(-2, 80))
```

```julia
using Distributions
using StatsPlots
using Random

dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
μ, σ = mean(dist), std(dist)
n = 320
L = 10^4
X = zeros(n)
MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
scatter(MV; ms=1, msw=0, alpha=0.5, label="",
    title="MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05]), n=$n",
    titlefontsize=12)
scatter!([(μ, σ^2)]; ms=3, label="")
plot!(; xlim=(-0.5, 4), ylim=(-2, 80))
```

```julia
using Distributions
using StatsPlots
using Random
using QuadGK
using LinearAlgebra
default(fmt=:png)

mystdmoment(dist::ContinuousUnivariateDistribution, k) = quadgk(x -> ((x - mean(dist))/std(dist))^k * pdf(dist, x), extrema(dist)...)[1]
myskewness(dist::MixtureModel{Univariate, Continuous}) = mystdmoment(dist, 3)
mykurtosis(dist::MixtureModel{Univariate, Continuous}) = mystdmoment(dist, 4) - 3
myskewness(dist::ContinuousUnivariateDistribution) = skewness(dist)
mykurtosis(dist::ContinuousUnivariateDistribution) = kurtosis(dist)

mydistname(dist) = replace(string(dist), r"\{.*\}"=>"")
function mydistname(dist::MixtureModel)
    c = components(dist)
    p = probs(dist)
    s = string(p[1]) * mydistname(c[1])
    for i in 2:length(p)
        s = s * "+" * string(p[i]) * mydistname(c[i])
    end
    s
end

function plot_samplemeanstd(dist, n; L=10^4, kwargs...)
    μ, σ = mean(dist), std(dist)
    sk, ku = myskewness(dist), mykurtosis(dist)
    
    X = rand(dist, n)
    MV = [(rand!(dist, X); (mean(X), var(X))) for _ in 1:L]
    
    A = [
        σ^2/n    sk/n*σ^3
        sk/n*σ^3 (ku/n + 2/(n-1))*σ^4
    ] |> Symmetric
    mvnormal = MvNormal([μ, σ^2], A)
    MVN = rand(mvnormal, L)
    MN, VN = MVN[1,:], MVN[2,:]
    
    xlim, ylim = extrema(MN), extrema(VN)
    
    P = scatter(MV; ms=1, msw=0, alpha=0.5, label="",
        title="$(mydistname(dist)), n=$n", titlefontsize=8,
        xlim, ylim)
    scatter!([(μ, σ^2)]; ms=3, label="")
    plot!(; kwargs...)
    
    Q = scatter(MVN[1,:], MVN[2,:]; ms=1, msw=0, alpha=0.5, label="",
        title="mvnormal approx by the central limit theorem", titlefontsize=8,
        xlim, ylim)
    scatter!([(μ, σ^2)]; ms=3, label="")
    plot!(; kwargs...)
    
    plot(P, Q; size=(800, 400))
end
```

```julia
dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
myskewness(dist), mykurtosis(dist)
```

```julia
dist = MixtureModel([Normal(), Normal(20,1)], [0.95, 0.05])
plot_samplemeanstd(dist, 20)
```

```julia
plot_samplemeanstd(dist, 40)
```

```julia
plot_samplemeanstd(dist, 80)
```

```julia
plot_samplemeanstd(dist, 160)
```

```julia
plot_samplemeanstd(dist, 320)
```

```julia

```
