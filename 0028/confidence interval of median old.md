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
    display_name: Julia 1.7.2
    language: julia
    name: julia-1.7
---

```julia
using Distributions
using StatsPlots
default(titlefontsize=12, fmt=:png)
using Roots
using Random
using StatsBase
using Memoization
```

```julia
name(dist::UnivariateDistribution) = replace(string(dist), r"{[^{.]*}"=>"")

function cdfordstat(dist, n, k, x)
    beta = Beta(k, n-k+1)
    cdf(beta, cdf(dist, x))
end

prevx(x::AbstractFloat) = prevfloat(x)
prevx(x::Integer) = x - 1
pmfordstat(dist, n, k, x) = cdfordstat(dist, n, k, x) - cdfordstat(dist, n, k, prevx(x))

function quantileordstat(dist, n, k, p)
    beta = Beta(k, n-k+1)
    quantile(dist, quantile(beta, p))
end

pmfmedian(dist, n, a) = pmfordstat(dist, n, (n+1)/2, a)
cdfmedian(dist, n, a) = cdfordstat(dist, n, (n+1)/2, a)
quantilemedian(dist, n, p) = quantileordstat(dist, n, (n+1)/2, p)
pvalmedian(dist, n, a) = min(1, 2cdfmedian(dist, n, a), 2(1 - cdfmedian(dist, n, prevx(a))))

function empiricaldist_allowrep(X) # very slow
    c = countmap(X)
    u = collect(keys(c))
    w = values(c) ./ length(X)
    DiscreteNonParametric(u, w) 
end
empiricaldist(X) = DiscreteNonParametric(X, fill(1/length(X), length(X)))

pmfmedian(X, a) = pmfmedian(empiricaldist(X), length(X), a)
cdfmedian(X, a) = cdfmedian(empiricaldist(X), length(X), a)
pvalmedian(X, a) = min(1, 2cdfmedian(X, a), 2(1 - cdfmedian(X, prevx(a))))

pmfmedian_allowrep(X, a) = pmfmedian(empiricaldist_allowrep(X), length(X), a)
cdfmedian_allowrep(X, a) = cdfmedian(empiricaldist_allowrep(X), length(X), a)
pvalmedian_allowrep(X, a) = min(1, 2cdfmedian_allowrep(X, a), 2(1 - cdfmedian(X, prevx(a))))

function cimedian_old(sampledist::DiscreteNonParametric, n; α = 0.05) # very slow
    f(x) = cdfmedian(sampledist, n, x)
    L = find_zeros(x -> f(x) - α/2,       extrema(sampledist)...)[begin]
    U = find_zeros(x -> f(x) - (1 - α/2), extrema(sampledist)...)[end]
    L, U
end

function cimedian(sampledist, n; α = 0.05)
    L = quantilemedian(sampledist, n, α/2)
    U = quantilemedian(sampledist, n, 1 - α/2)
    L, U
end

function cimedian(X; α = 0.05)
    n = length(X)
    sampledist = empiricaldist(X)
    cimedian(sampledist, n; α)
end

function cimedian_allowrep(X; α = 0.05)
    n = length(X)
    sampledist = empiricaldist_allowrep(X)
    cimedian(sampledist, n; α)
end

function plot_randmedianci(dist, n; α = 0.05)
    X = rand(dist, n)
    sampledist = empiricaldist(X)
    a, b = minmax(median(dist), median(X))
    s = max(std(dist), std(X))/√n
    xlim = (min(a,b) - 5s, max(a,b) + 5s)
    ci = cimedian(X; α) |> collect
    plot(; legend=:topleft)
    plot!(x -> cdfmedian(dist, n, x), xlim...; label="true")
    plot!(x -> cdfmedian(sampledist, n, x), xlim...; label="bootstrap estimation")
    vline!([median(dist)]; label="true median", c=1, ls=:dash)
    vline!([median(X)]; label="sample median", c=2, ls=:dash)
    plot!(ci, fill(α, 2); label="confidence interval", c=:red, alpha=0.7, lw=3)
    plot!(; ytick=0:0.05:1)
    title!("$(name(dist)), n=$n")    
end
```

```julia
x ⪅ y = x < y || x ≈ y
x ⪉ y = x < y && !(x ≈ y)

function myquantile(d::DiscreteNonParametric, q::Real)
    0 <= q <= 1 || throw(DomainError())
    x = support(d)
    p = probs(d)
    k = length(x)
    i = 1
    cp = p[1]
    while cp ⪉ q && i < k #Note: is i < k necessary?
        i += 1
        @inbounds cp += p[i]
    end
    i == k && return x[i]
    cp ⪅ q && return (p[i]*x[i] + p[i+1]*x[i+1])/(p[i] + p[i+1])
    i > 1 && return (p[i-1]*x[i-1] + p[i]*x[i])/(p[i-1] + p[i])
    x[i]
end

mymedian(d::DiscreteNonParametric) = myquantile(d, 0.5)
```

```julia
X = rand(Normal(2, 3), 10)
@show median(X)
@show empiricaldist(X) |> median
@show empiricaldist(X) |> mymedian
@show sort(X)[5:6]
@show cdf.(empiricaldist(X), sort(X)[5:6]);
```

```julia
X = rand(Normal(2, 3), 11)
@show median(X)
@show empiricaldist(X) |> median
@show empiricaldist(X) |> mymedian
@show sort(X)[5:7]
@show cdf.(empiricaldist(X), sort(X)[5:7]);
```

```julia
plot_randmedianci(Normal(2, 3), 10; α = 0.05)
```

```julia
plot_randmedianci(Normal(2, 3), 20; α = 0.05)
```

```julia
plot_randmedianci(Normal(2, 3), 40; α = 0.05)
```

```julia
plot_randmedianci(Normal(2, 3), 80; α = 0.05)
```

```julia
plot_randmedianci(Gamma(2, 3), 10; α = 0.05)
```

```julia
plot_randmedianci(Gamma(2, 3), 20; α = 0.05)
```

```julia
plot_randmedianci(Gamma(2, 3), 40; α = 0.05)
```

```julia
plot_randmedianci(Gamma(2, 3), 80; α = 0.05)
```

```julia
plot_randmedianci(Exponential(), 10; α = 0.05)
```

```julia
plot_randmedianci(Exponential(), 20; α = 0.05)
```

```julia
plot_randmedianci(Exponential(), 40; α = 0.05)
```

```julia
plot_randmedianci(Exponential(), 80; α = 0.05)
```

```julia
dist, n = Normal(2, 3), 100
X = rand(dist, n)
plot(x -> pvalmedian(X, x), -1, 5; label="pval")
```

```julia
dist, n = Gamma(2, 3), 100
X = rand(dist, n)
plot(x -> pvalmedian(X, x), 3, 8; label="pval")
```

```julia
dist, n = Exponential(), 100
X = rand(dist, n)
plot(x -> pvalmedian(X, x), 0, 1.5; label="pval")
```

```julia
function sim_mediantest(dist, n; L = 10^5)
    a = median(dist)
    pval = Vector{Float64}(undef, L)
    tmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, tmp[Threads.threadid()])
        pval[i] = pvalmedian(X, a)
    end
    pval
end

function plot_mediantest(dist, n; L = 10^5)
    pval = sim_mediantest(dist, n; L)
    m = median(dist)
    s = std(dist)/√n
    plot(; legend=false)
    plot!(a -> ecdf(pval)(a), 0, 0.1; label="ecdf of pvalues")
    plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot)
    plot!(; xtick=0:0.01:1, ytick=0:0.01:1)
    plot!(; size=(400, 400))
end
```

```julia
n = 10
P1 = plot_mediantest(Normal(2, 3), n)
P2 = plot_mediantest(Gamma(2, 3), n)
P3 = plot_mediantest(LogNormal(), n)
plot(P1, P2, P3; size=(900, 300), layout=(1, 3), tickfontsize=6)
```

```julia
n = 20
P1 = plot_mediantest(Normal(2, 3), n)
P2 = plot_mediantest(Gamma(2, 3), n)
P3 = plot_mediantest(LogNormal(), n)
plot(P1, P2, P3; size=(900, 300), layout=(1, 3), tickfontsize=6)
```

```julia
n = 30
P1 = plot_mediantest(Normal(2, 3), n)
P2 = plot_mediantest(Gamma(2, 3), n)
P3 = plot_mediantest(LogNormal(), n)
plot(P1, P2, P3; size=(900, 300), layout=(1, 3), tickfontsize=6)
```

```julia
n = 40
P1 = plot_mediantest(Normal(2, 3), n)
P2 = plot_mediantest(Gamma(2, 3), n)
P3 = plot_mediantest(LogNormal(), n)
plot(P1, P2, P3; size=(900, 300), layout=(1, 3), tickfontsize=6)
```

```julia
n = 80
P1 = plot_mediantest(Normal(2, 3), n)
P2 = plot_mediantest(Gamma(2, 3), n)
P3 = plot_mediantest(LogNormal(), n)
plot(P1, P2, P3; size=(900, 300), layout=(1, 3), tickfontsize=6)
```

```julia
n = 160
P1 = plot_mediantest(Normal(2, 3), n)
P2 = plot_mediantest(Gamma(2, 3), n)
P3 = plot_mediantest(LogNormal(), n)
plot(P1, P2, P3; size=(900, 300), layout=(1, 3), tickfontsize=6)
```

観察: 第一種の過誤が起こる確率は母集団分布によらずに決まっている.

```julia
dist, n = Gamma(2, 3), 10
X = rand(dist, n)
Y = -4.0:n-5
@show cdfX = cdfmedian.(Ref(X), sort(X))
@show cdfY = cdfmedian.(Ref(Y), Y)
@show cdfX == cdfY
println()
@show pmfX = pmfmedian.(Ref(X), sort(X))
@show sum(pmfX)
@show pmfY = pmfmedian.(Ref(Y), Y)
@show sum(pmfY)
@show pmfX == pmfY
println()
@show pvalX = pvalmedian.(Ref(X), sort(X))
@show pvalY = pvalmedian.(Ref(Y), Y)
@show pvalX == pvalY;
```

このように, 任意のサイズ $n$ のサンプル $X_1<X_2<\cdots<X_n$ から得られるP値函数の $X_i$ における値と, $1<2<\cdots<n$ から得られるP値函数の $i$ における値は等しくなる.

```julia
X = Float64[5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 2, 1, 1, 1, 1]
@show X
@show Xdist = empiricaldist_allowrep(X)
@show n = 10
P1 = plot(x -> cdfmedian(Xdist, n, x), 0, 6; label="", ytick=0:0.1:1)
P2 = plot(x -> pvalmedian(Xdist, n, x), 0, 6; label="", ytick=0:0.1:1)
plot(P1, P2; size=(800, 300))
```

```julia
X = Float64[5, 5, 5, 5, 5, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 2, 1, 1, 1, 1]
Xdist = Xdist = empiricaldist_allowrep(X)
n = 15
@show X
@show Xdist
@show n
P1 = plot(x -> cdfmedian(Xdist, n, x), 0, 6; label="", ytick=0:0.1:1)
P2 = plot(x -> pvalmedian(Xdist, n, x), 0, 6; label="", ytick=0:0.1:1)
plot(P1, P2; size=(800, 300))
```

## RのDescToolsによる中央値の信頼区間の計算との比較

https://stats.stackexchange.com/questions/502977/confidence-interval-for-median-which-is-more-appropriate-bootstrap-or-binom-ex

![iVnfe.png](attachment:956643d6-7cb2-482e-a0dc-96efe43a57cb.png)

```julia
using RCall
R"library(DescTools)"
```

```julia
X = data1 = [8,    7,  8,  9.5,    1,  20, 8,  7.5,    3,  20.5,   2.5,    5.5,    15.5,   2,  4,  1,
    17,   2,  3.5,    8.5,    8.5,    2.5,    11, 4,  10.5,   7.5,    12, 5,  16.5,   8.5]
Y = X + [√eps()*randn() for _ in 1:length(X)]
@show cimedian_allowrep(X)
@show cimedian(Y)
println()
@rput data1
@show R"""MedianCI(data1, na.rm = TRUE, method = "boot")"""
@show R"""MedianCI(data1, na.rm = TRUE, method = "exact")""";
```

```julia
X = data2 = [7.1,  32.0,   3.8,    1.6,    19.6,   6.0,    7.2,    14.9,   0,  2.0,    5.7,
          19.4, 13.1,   15.5,   11.3,   9.6,    13.9,   5.6,    12.6,   1.0,    1.9,
          8.1,  15.9,   0.8,    6.1,    8.1,    18.0,   4.6,    5.5,    15.6]
Y = X + [√eps()*randn() for _ in 1:length(X)]
@show cimedian_allowrep(X)
@show cimedian(Y)
println()
@rput data2
@show R"""MedianCI(data2, na.rm = TRUE, method = "boot")"""
@show R"""MedianCI(data2, na.rm = TRUE, method = "exact")""";
```

```julia
X = data3 = [16.1, 10.4,   0.5,    12.2,   7.2,    1.7,    21.6,   6.3,    0.8,    3.2,    12.6,   20.0,   3.4, 7.3,   3.5,
          7.5,  15.8, 4.7, 8.3, 11.9,   1.6,    9.0, 8.6,   11.7,   8.1, 5.8, 3.3,  7.9,    7.0,    8.5]
Y = X + [√eps()*randn() for _ in 1:length(X)]
@show cimedian_allowrep(X)
@show cimedian(Y)
println()
@rput data3
@show R"""MedianCI(data3, na.rm = TRUE, method = "boot")"""
@show R"""MedianCI(data3, na.rm = TRUE, method = "exact")""";
```

https://github.com/cran/DescTools/blob/d5e096ee9abf4640703dfba45d1ed56b5ab10253/R/StatsAndCIs.r#L3690

* https://www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf
* https://www.scribd.com/presentation/75941305/Confidence-Interval-for-Median-Based-on-Sign-Test
* https://stat.ethz.ch/pipermail/r-help/2003-September/039636.html

このノートで実装した中央値の信頼区間は "SAS-way" と同じもののようだ.

```julia
p = [1, 2, 3, 0, 4, 5, 7, 2, 1, 2, 3]
p = p/sum(p)
mix = MixtureModel([Uniform(k, k+1) for k in 1:length(p)], p)
```

```julia
quantile(mix, 0.2)
```

```julia
cimedian(mix, 30)
```

```julia
n = 30
P1 = plot(x -> cdfmedian(mix, n, x), 5, 8.5; label="", ytick=0:0.1:1)
P2 = plot(x -> pvalmedian(mix, n, x), 5, 8.5; label="", ytick=0:0.1:1)
plot(P1, P2; size=(800, 300))
```

```julia

```
