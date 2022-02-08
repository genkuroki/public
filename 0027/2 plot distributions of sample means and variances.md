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

$
\newcommand\var{\operatorname{var}}
\newcommand\cov{\operatorname{cov}}
$


__記号の定義__

Let $X_1,\ldots,X_n$ be a sample with size $n$ of a distribution with mean $\mu$ and variance $\sigma^2$ and set

$$
\begin{aligned}
&
\bar{X} = \frac{1}{n}\sum_{i=1}^n X_i,
\\ &
S^2 = \frac{1}{n-1}\sum_{i=1}^n (X_i - \bar{X})^2,
\\ &
Z = \frac{\sqrt{n}(X - \mu)}{\sigma}
\quad (\text{expected to follow the standard normal distribution}),
\\ &
W = \sqrt{n}\frac{S^2}{\sigma^2},
\\ &
X^2 = \frac{(n-1)S^2}{\sigma^2} = \frac{1}{\sigma^2}\sum_{i=1}^n (X_i - \bar{X})^2
\quad (\text{expected to follow the chi-sauared distribution with df = $n-1$}),
\\ &
T = \frac{\sqrt{n}(X - \mu)}{S} = \frac{\sigma}{S}Z
\quad (\text{expected to follow the t-distribution with df = $n-1$}).
\end{aligned}
$$

Then $T^2$ is expected to follow the F-distribution with dfs = ($1$, $n-1$).


$Z = \sqrt{n}(\bar{X} - \mu)/\sigma$ and $W = \sqrt{n}\;S^2/\sigma^2$ satisfy
$$
\begin{aligned}
\\ &
E[Z] = 0, \quad E[W] = \sqrt{n},
\\ &
\var(Z) = 1,
\quad
\cov(Z, W) = \kappa_3
\quad
\var(W) = \kappa_4 + 2 + \frac{2}{n-1}
\\ &
E[Z^2] = 1, \quad E[W] = \sqrt{n},
\\ &
\var(Z^2) = \frac{\kappa_4}{n} + 2,
\quad
\cov(Z^2, W) = \frac{\kappa_4}{\sqrt{n}},
\quad
\var(W) = \kappa_4 + 2 + \frac{2}{n-1}.
\end{aligned}
$$


$Z = \sqrt{n}(\bar{X} - \mu)/\sigma$ and $R^2 = S^2/\sigma^2$ satisfy
$$
\begin{aligned}
\\ &
E[Z] = 0, \quad E[R^2] = 1,
\\ &
\var(Z) = 1,
\quad
\cov(Z, R^2) = \frac{\kappa_3}{\sqrt{n}}
\quad
\var(R^2) = \frac{\kappa_4}{n} + \frac{2}{n-1}
\\ &
E[Z^2] = 1, \quad E[R^2] = 1,
\\ &
\var(Z^2) = \frac{\kappa_4}{n} + 2,
\quad
\cov(Z^2, R^2) = \frac{\kappa_4}{n},
\quad
\var(R^2) = \frac{\kappa_4}{n} + \frac{2}{n-1}.
\end{aligned}
$$


__この Part 2 でのプロット__

$(Z^2, S^2/\sigma^2)$ をプロットしたい.


__Part 1 でのプロット__

次の４つをプロットした:

1. $(Z, X^2)$ のサンプルの散布図
2. $X^2$ の値を決めたときの $Z$ の条件付き確率分布の密度函数 $p(z|x^2)$ のサンプルからの推定結果
3. $X^2$ の分布と自由度 $n-1$ のカイ二乗分布の比較
4. $T^2$ の分布と自由度 $(1, n-1)$ のF分布の比較

$T$ そのものの分布ではなく、その二乗の分布を見ることに注意。検定で使われるのは $T$ の絶対値だけなのであった。 

2におけるシアンの線はZの値の下位５％と上位95%の部分の境界である。推定値なので少しギザギザしてしまっている。

$Z$ と $X^2$ が確率変数として独立な場合にはシアンの線が縦軸に平行になる。

```julia
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=8, legendfontsize=7, guidefontsize=8, tickfontsize=6)
using Random
using KernelDensity
using QuadGK
using StatsBase

safediv(x, y) = y == 0 ? y : x/y

function myquantile(y, Y, Z, p; δ=log2(y)/4)
    mask = @. y - δ < Y < y + δ
    !any(mask) && return NaN
    quantile(Z[mask], p)
end

function mymean(y, Y, Z, p; δ=log2(y)/4)
    mask = @. y - δ < Y < y + δ
    !any(mask) && return NaN
    mean(Z[mask])
end

myskewness(dist) = skewness(dist)
myskewness(dist::MixtureModel) = _myskewness(dist)
function _myskewness(dist)
    μ, σ = mean(dist), std(dist)
    f(x) = ((x - μ)/σ)^3 * pdf(dist, x)
    quadgk(f, extrema(dist)...)[1]
end

mykurtosis(dist) = kurtosis(dist)
mykurtosis(dist::MixtureModel) = _mykurtosis(dist)
function _mykurtosis(dist)
    μ, σ = mean(dist), std(dist)
    f(x) = ((x - μ)/σ)^4 * pdf(dist, x)
    quadgk(f, extrema(dist)...)[1] - 3
end

rd(x; digits=4) = round(x; digits)

name(dist) = replace(string(dist), r"\{.*\}"=>"")
name(dist::MixtureModel) = _name(dist)
function _name(dist::MixtureModel)
    c = components(dist)
    p = probs(dist)
    s = string(p[1]) * name(c[1])
    for i in 2:length(p)
        s = s * "+" * string(p[i]) * name(c[i])
    end
    s
end
```

```julia
function mcsim(;
        dist = Beta(0.2, 0.3),
        n = 20,
        L = 10^6,
    )
    μ, σ = mean(dist), std(dist)
    Z = Vector{Float64}(undef, L) # expected to follow Normal(1,0)
    X² = similar(Z)               # expected to follow Chisq(n-1)
    tmp = [Vector{Float64}(undef, n) for i in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, tmp[Threads.threadid()])
        X̄ = mean(X) # sample mean
        S² = var(X) # sample unbiased variance
        Z[i] = √n * (X̄ - μ) / σ
        X²[i] = (n - 1) * S²/σ^2
    end
    T = @. Z / √(X²/(n - 1)) # expected to follow TDist(n-1)
    T² = @. T .^ 2           # expected to follow FDist(1, n-1)
    (; dist, n, Z, X², T, T²)
end

function plot_samplestats(;
        dist = Normal(1, 2),
        n = 10,
        L = 10^6,
        scattermax = 10^5,
        T²span = (0.75, 0.999),
        Zspan = (0.001, 0.999), 
        X²span = (0.001, 0.999), 
        kwargs...
    )
    (; dist, n, Z, X², T, T²) = mcsim(; dist, n, L)
    distname = name(dist)
    sk = myskewness(dist)
    ku = mykurtosis(dist)
    fdist = FDist(1, n-1)
    
    println(L, " samples with size n = ", n, " of ", distname)
    println("skewness, kurtosis = ", rd(sk), ", ", rd(ku))
    for p in (0.95,)# 0.99)
        println("√quantile(T², $p) = ", rd(√quantile(T², p)), ",  ")
        println("P(|T| > √quantile(FDist(1, $(n-1)), $p) = ", rd(√quantile(fdist, p)), ") = ",  rd(1 - ecdf(T²)(quantile(fdist, p))))
    end
    
    Zlim = quantile.(Ref(Z), Zspan)
    X²lim = quantile.(Ref(X²), X²span)
    
    kdeX² = InterpKDE(kde(X²))
    kdeZX² = InterpKDE(kde((Z, X²)))
    g(z, x²) = safediv(pdf(kdeZX², z, x²), pdf(kdeX², x²))
    kdeT² = InterpKDE(kde(T²))
    h(t²) = pdf(kdeT², t²)
    
    P = plot(; colorbar=false)
    plot!(; xlabel="Z = √n(X̄ - μ)/σ", ylabel="X² = (n - 1)S²/σ²", xlim=Zlim, ylim=X²lim)
    scatter!(Z[1:min(end, scattermax)], X²[1:min(end, scattermax)]; alpha=0.3, msw=0, ms=1, label="")
    vline!([0]; label="", ls=:dot, c=:red)
    hline!([n-1]; label="", ls=:dot, c=:red)
    title!("$distname, n=$n")
    
    z = range(Zlim..., 200)
    x² = range(X²lim..., 200)
    z05 = (x² -> myquantile(x², X², Z, 0.05)).(x²)
    z95 = (x² -> myquantile(x², X², Z, 0.95)).(x²)
    Q = plot(; colorbar=false)
    plot!(; xlabel="Z = √n(X̄ - μ)/σ conditioned by X²", ylabel="X² = (n - 1)S²/σ²", xlim=Zlim, ylim=X²lim)
    heatmap!(z, x², g)
    plot!(z05, x²; label="", c=:cyan)
    plot!(z95, x²; label="", c=:cyan)
    vline!([0]; label="", ls=:dot, c=:pink)
    hline!([n-1]; label="", ls=:dot, c=:pink)
    title!("p(z|x²)")
    
    chisqdist = Chisq(n-1)
    xlim = (max(0.0, quantile(X², 0.005) - 10), max(quantile(X², 0.995), quantile(chisqdist, 0.999)))
    R = plot(; xlabel="X² = (n - 1)S²/σ²", ylabel="density", xlim)
    histogram!(X²; norm=true, alpha=0.3, bin=range(xlim..., 100), label="X²")
    plot!(chisqdist, xlim...; label="Chisq($(n-1))", lw=1.5)
    vline!([n-1]; label="", ls=:dot, c=:black)
    title!("sample of X² = (n - 1)S²/σ²")
    
    xlim = quantile.(Ref(T²), T²span)
    bin = range(0, last(xlim), round(Int, 100last(xlim)/(last(xlim) - first(xlim))))
    ymax = maximum(x -> max(pdf(fdist, x), h(x)), range(xlim..., 100))
    ylim = (-0.03ymax, 1.05ymax)
    S = plot(; xlabel="T² where T = √n(X̄ - μ)/S", ylabel="density", xlim, ylim)
    histogram!(T²; norm=true, alpha=0.3, bin=bin, label="T²")
    plot!(fdist, xlim...; label="FDist(1, $(n-1))", lw=1.5)
    vline!([quantile(T², 0.95)]; label="95% line of T²", c=1, ls=:dot)
    vline!([quantile(fdist, 0.95)]; label="95% line of FDist", c=2, ls=:dot)
    title!("tail (> $(100first(T²span))%) of sample of T²")
    
    plot(P, Q, R, S; size=(800, 600))
    plot!(leftmargin=3Plots.mm, bottommargin=3Plots.mm, kwargs...)
end
```

```julia
function simZR(;
        dist = Uniform(),
        n = 20,
        L = 10^6,
    )
    μ, σ = mean(dist), std(dist)
    Z = Vector{Float64}(undef, L) # Z² is expected to follow Chisq(1)/1.        mean(Z²) ≈ 1, var(Z²) ≈ 2
    R² = similar(Z)               # R² is expected to follow Chisq(n-1)/(n-1).  mean(R²) ≈ 1, var(R²) ≈ 2/(n-1)
    tmp = [Vector{Float64}(undef, n) for i in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, tmp[Threads.threadid()])
        X̄ = mean(X) # sample mean
        S² = var(X) # sample unbiased variance
        Z[i] = √n * (X̄ - μ) / σ
        R²[i] = S² / σ^2
    end
    Z² = @. Z^2
    T² = @. Z² / R² # expected to follow FDist(1, n-1)
    (; dist, n, Z, Z², R², T²)
end

function plotZR(;
        dist = Uniform(),
        n = 10,
        L = 10^6,
        scattermax = 10^5,
        Zspan = (0.001, 0.999), 
        Z²span = (0.001, 0.999), 
        R²span = (0.001, 0.999), 
        T²span = (0.001, 0.999), 
        kwargs...
    )
    (; dist, n, Z, Z², R², T²) = simZR(; dist, n, L)
    distname = name(dist)
    κ₃ = myskewness(dist)
    κ₄ = mykurtosis(dist)
    fdist = FDist(1, n-1)
   
    println(L, " samples with size n = ", n, " of ", distname)
    println("skewness, kurtosis = ", rd(κ₃), ", ", rd(κ₄))
    for p in (0.95,)# 0.99)
        println("√quantile(T², $p) = ", rd(√quantile(T², p)), ",  ")
        println("P(|T| > √quantile(FDist(1, $(n-1)), $p) = ", rd(√quantile(fdist, p)), ") = ",  rd(1 - ecdf(T²)(quantile(fdist, p))))
    end
     
    @show Z²lim = quantile.(Ref(Z²), Z²span)
    @show R²lim = quantile.(Ref(R²), R²span)

    @show mean(Z²), std(Z²), √(κ₄/n + 2)
    @show mean(R²), std(R²), √(κ₄/n + 2/(n-1))
    @show cov(Z², R²), κ₄/n

    kdeR² = InterpKDE(kde(R²))
    kdeZ²R² = InterpKDE(kde((Z², R²)))
    g(z², r²) = √(max(eps(), pdf(kdeZ²R², z², r²)) / max(eps(), pdf(kdeR², r²)))

    kdeT²R² = InterpKDE(kde((T², R²)))
    h(z², r²) = √(max(eps(), pdf(kdeT²R², z², r²)) / max(eps(), pdf(kdeR², r²)))
    #h(z², r²) = max(eps(), pdf(kdeT²R², z², r²)) / max(eps(), pdf(kdeR², r²))

    kdeT² = InterpKDE(kde(T²))
    h(t²) = pdf(kdeT², t²)
    
    P1 = plot(; colorbar=false)
    plot!(; xlabel="Z² = n(X̄ - μ)²/σ²", ylabel="R² = S²/σ²", xlim=Z²lim, ylim=R²lim)
    scatter!(Z²[1:min(end, scattermax)], R²[1:min(end, scattermax)]; alpha=0.3, msw=0, ms=1, label="")
    vline!([1]; label="", ls=:dot, c=:red)
    hline!([1]; label="", ls=:dot, c=:red)
    title!("$distname, n=$n")
    
    z² = range(Z²lim..., 200)
    t² = range(Z²lim..., 200)
    r² = range(R²lim..., 200)
    δ = (last(R²lim) - first(R²lim))/200
    
    λ₁ = cov(Z², R²)/var(R²)
    @show λ₁
    z²meanexpected = (r² -> r² ^ λ₁).(r²)
    z²mean = (r² -> mymean(r², R², Z², 0.95; δ)).(r²)
    z²95 = (r² -> myquantile(r², R², Z², 0.95; δ)).(r²)
   
    P2 = plot(; colorbar=false, legend=:right)
    plot!(; xlabel="Z² = n(X̄ - μ)²/σ² conditioned by R²", ylabel="R² = S²/σ²", xlim=Z²lim, ylim=R²lim)
    heatmap!(z², r², g)
    plot!(z²meanexpected, r²; label="r² ^ λ₁", c=:yellow)
    plot!(z²mean, r²; label="mean", c=:red, ls=:dash)
    plot!(z²95, r²; label="95%", c=:cyan)
    vline!([1]; label="", ls=:dot, c=:pink)
    hline!([1]; label="", ls=:dot, c=:pink)
    title!("√p(z²|r²)")
    
    λ₂ = √((κ₄/n + 2/(n-1)) / (2/(n-1)))
    @show -1 + λ₁
    @show - λ₂
    t²meanexpected1 = (r² -> r² ^ (- 1 + λ₁)).(r²)
    t²meanexpected2 = (r² -> r² ^ (- λ₂)).(r²)
    t²mean = (r² -> mymean(r², R², T², 0.50; δ)).(r²)
    t²95 = (r² -> myquantile(r², R², T², 0.95; δ)).(r²)
 
    P3 = plot(; colorbar=false, legend=:right)
    plot!(; xlabel="T² = n(X̄ - μ)²/S² conditioned by R²", ylabel="R² = S²/σ²", xlim=Z²lim, ylim=R²lim)
    #heatmap!(t², r², h)
    plot!(t²mean, r²; label="mean", c=:red)
    #plot!(t²meanexpected1, r²; label="", c=:yellow, ls=:dash)
    plot!(t²meanexpected2, r²; label="", c=:orange, ls=:dashdot)
    plot!(t²95, r²; label="95%", c=:cyan)
    hline!([1]; label="", ls=:dot, c=:pink)
    title!("√p(t²|r²)")
    
    P4 = plot(t², h)
    plot!(fdist, 0.1, extrema(t²)[2])
    
#     chisqdist = Chisq(n-1)
#     xlim = (max(0.0, quantile(X², 0.005) - 10), max(quantile(X², 0.995), quantile(chisqdist, 0.999)))
#     R = plot(; xlabel="X² = (n - 1)S²/σ²", ylabel="density", xlim)
#     histogram!(X²; norm=true, alpha=0.3, bin=range(xlim..., 100), label="X²")
#     plot!(chisqdist, xlim...; label="Chisq($(n-1))", lw=1.5)
#     vline!([n-1]; label="", ls=:dot, c=:black)
#     title!("sample of X² = (n - 1)S²/σ²")
    
#     xlim = quantile.(Ref(T²), T²span)
#     bin = range(0, last(xlim), round(Int, 100last(xlim)/(last(xlim) - first(xlim))))
#     ymax = maximum(x -> max(pdf(fdist, x), h(x)), range(xlim..., 100))
#     ylim = (-0.03ymax, 1.05ymax)
#     S = plot(; xlabel="T² where T = √n(X̄ - μ)/S", ylabel="density", xlim, ylim)
#     histogram!(T²; norm=true, alpha=0.3, bin=bin, label="T²")
#     plot!(fdist, xlim...; label="FDist(1, $(n-1))", lw=1.5)
#     vline!([quantile(T², 0.95)]; label="95% line of T²", c=1, ls=:dot)
#     vline!([quantile(fdist, 0.95)]; label="95% line of FDist", c=2, ls=:dot)
#     title!("tail (> $(100first(T²span))%) of sample of T²")
    
    plot(P1, P2, P3, P4; size=(800, 600))
    plot!(leftmargin=3Plots.mm, bottommargin=3Plots.mm, kwargs...)
end
```

```julia
plotZR(; dist = Normal(1,2), n = 10)
```

```julia
plotZR(; dist = Normal(1,2), n = 40)
```

```julia
plotZR(; dist = Uniform(), n = 10)
```

```julia
plotZR(; dist = Uniform(), n = 40)
```

```julia
plotZR(; dist = Uniform(), n = 160)
```

```julia
plotZR(; dist = Uniform(), n = 640)
```

```julia
plotZR(; dist=Beta(0.2, 0.2), n = 10)
```

```julia
plotZR(; dist=Beta(0.2, 0.2), n = 20)
```

```julia
plotZR(; dist=Beta(0.2, 0.2), n = 40)
```

```julia
plotZR(; dist=Beta(0.2, 0.2), n = 160)
```

```julia
plotZR(; dist=TDist(5), n = 10)
```

```julia
plotZR(; dist=TDist(5), n = 40)
```

```julia
plotZR(; dist=TDist(5), n = 160)
```

```julia
plotZR(; dist=TDist(5), n = 640)
```

```julia
plotZR(; dist=Laplace(), n = 10)
```

```julia
plotZR(; dist=Laplace(), n = 40)
```

```julia
plotZR(; dist=Laplace(), n = 160)
```

```julia
plotZR(; dist=Laplace(), n = 640)
```

```julia

```

```julia

```

```julia

```

```julia

```

```julia

```

```julia

```

```julia

```

```julia
plotZR(; dist=Gamma(5,1), n = 10)
```

```julia
plotZR(; dist=Gamma(5,1), n = 40)
```

```julia
plotZR(; dist=Gamma(5,1), n = 160)
```

```julia
plotZR(; dist=Gamma(5,1), n = 640)
```

```julia
plotZR(; dist=Gamma(20,1), n = 10)
```

```julia
plotZR(; dist=Gamma(20,1), n = 40)
```

```julia
plotZR(; dist=Gamma(20,1), n = 160)
```

```julia
plotZR(; dist=Gamma(20,1), n = 640)
```

```julia
plotZR(; dist=Exponential(), n = 10)
```

```julia
plotZR(; dist=Exponential(), n = 40)
```

```julia
plotZR(; dist=Exponential(), n = 160)
```

```julia
plotZR(; dist=Exponential(), n = 640)
```

```julia

```

```julia

```

```julia

```

```julia

```

```julia

```

```julia
dist = Normal(1, 2); plot(dist, -7, 9; label="", title=name(dist), titlefontsize=12)
```

```julia
plot_samplestats(dist = Normal(1, 2), n = 10)
```

```julia
dist = Uniform(); plot(dist, -0.2, 1.2; label="", title=name(dist), titlefontsize=12)
```

```julia
plot_samplestats(dist = Uniform(), n = 10)
```

```julia
plot_samplestats(dist = Uniform(), n = 20)
```

```julia
plot_samplestats(dist = Uniform(), n = 40)
```

```julia
plot_samplestats(dist = Uniform(), n = 80)
```

```julia
plot_samplestats(dist = Uniform(), n = 160)
```

```julia
plot_samplestats(dist = Uniform(), n = 320)
```

```julia
plot_samplestats(dist = Uniform(), n = 640)
```

```julia
plot_samplestats(dist = Uniform(), n = 1280)
```

```julia
plot_samplestats(dist = Uniform(), n = 2560)
```

```julia
dist = Beta(0.2, 0.2); plot(dist, 0.005, 0.995; label="", title=name(dist), titlefontsize=12)
```

```julia
plot_samplestats(dist = Beta(0.2, 0.2), n = 10)
```

```julia
plot_samplestats(dist = Beta(0.2, 0.2), n = 20)
```

```julia
dist = Laplace(); plot(dist, -6, 6; label=name(dist))
plot!(Normal(), -6, 6; label="Normal()", ls=:dash)
```

```julia
plot_samplestats(dist = Laplace(), L=10^6)
```

```julia
dist = TDist(4.01); plot(dist, -6, 6; label=name(dist))
plot!(Normal(), -6, 6; label="Normal()", ls=:dash)
```

```julia
plot_samplestats(dist = TDist(4.01), L=10^6, X²span = (0.025, 0.975))
```

```julia
plot_samplestats(dist = TDist(4.01), n = 20, L=10^6, X²span = (0.025, 0.975))
```

```julia
plot_samplestats(dist = TDist(4.01), L=10^6, X²span = (0.025, 0.975))
```

```julia
plot_samplestats(dist = TDist(4.01), L=10^6, X²span = (0.025, 0.975))
```

```julia
plot_samplestats(dist = TDist(4.01), n = 160, L=10^6, X²span = (0.025, 0.975))
```

```julia
plot_samplestats(dist = TDist(4.01), n = 320, L = 10^6, X²span = (0.025, 0.975))
```

```julia
dist = TDist(5); plot(dist, -6, 6; label=name(dist))
plot!(Normal(), -6, 6; label="Normal()", ls=:dash)
```

```julia
plot_samplestats(dist = TDist(5), n = 10, L = 10^6, X²span = (0.01, 0.99))
```

```julia
dist = TDist(10); plot(dist, -6, 6; label=name(dist))
plot!(Normal(), -6, 6; label="Normal()", ls=:dash)
```

```julia
plot_samplestats(dist = TDist(10), n = 10)
```

```julia
plot_samplestats(dist = TDist(10), n = 20)
```

```julia
dist = MixtureModel([Normal(), Normal(10, 1)], [0.5, 0.5]); plot(dist, -6, 16; label="", title=name(dist))
```

```julia
plot_samplestats(; dist, n = 10, X²span = (0.01, 0.999), T²span = (0.50, 0.99))
```

```julia
plot_samplestats(; dist, n = 20, X²span = (0.01, 0.999), T²span = (0.50, 0.99))
```

```julia
dist = MixtureModel([Laplace(), Laplace(20, 1)], [0.5, 0.5]); plot(dist, -8, 28; label="", title=name(dist))
```

```julia
plot_samplestats(; dist, n = 10, X²span = (0.01, 0.999), T²span = (0.50, 0.99))
```

```julia
plot_samplestats(; dist, n = 20, X²span = (0.01, 0.999), T²span = (0.50, 0.99))
```

```julia
dist = Exponential(); plot(dist, 0, 6; label="", title=name(dist), titlefontsize=12)
```

```julia
plot_samplestats(dist = Exponential(), n = 10)
```

```julia
plot_samplestats(dist = Exponential(), n = 20)
```

```julia
plot_samplestats(dist = Exponential(), n = 40)
```

```julia
plot_samplestats(dist = Exponential(), n = 80)
```

```julia
plot_samplestats(dist = Exponential(), n = 160)
```

```julia
dist = Gamma(5, 1); plot(dist, 0, 15; label="", title=name(dist), titlefontsize=12)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 10)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 20)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 40)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 80)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 160)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 320)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 640)
```

```julia
plot_samplestats(dist = Gamma(5, 1), n = 1280)
```

```julia
dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05])
plot(dist, -4, 14; label="", title=name(dist), titlefontsize=12)
```

```julia
plot_samplestats(dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05]), n = 10)
```

```julia
plot_samplestats(dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05]), n = 20)
```

```julia
plot_samplestats(dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05]), n = 40)
```

```julia
plot_samplestats(dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05]), n = 80)
```

```julia
plot_samplestats(dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05]), n = 160)
```

```julia
plot_samplestats(dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05]), n = 320)
```

```julia
plot_samplestats(dist = MixtureModel([Normal(), Normal(10,1)], [0.95, 0.05]), n = 640)
```

```julia
dist = LogNormal()
plot(dist, 0, 10; label="", title=name(dist), titlefontsize=12)
```

```julia
plot_samplestats(dist = LogNormal(), n = 10)
```

```julia
plot_samplestats(dist = LogNormal(), n = 20)
```

```julia
plot_samplestats(dist = LogNormal(), n = 40)
```

```julia
plot_samplestats(dist = LogNormal(), n = 80)
```

```julia
plot_samplestats(dist = LogNormal(), n = 160)
```

```julia
plot_samplestats(dist = LogNormal(), n = 320)
```

```julia
plot_samplestats(dist = LogNormal(), n = 640)
```

```julia
plot_samplestats(dist = LogNormal(), n = 1280)
```

```julia
plot_samplestats(dist = LogNormal(), n = 2560)
```

```julia
plot_samplestats(dist = LogNormal(), n = 5120)
```

```julia
@time plot_samplestats(dist = LogNormal(), n = 10240)
```

```julia
@time plot_samplestats(dist = LogNormal(), n = 20480)
```

```julia
@time plot_samplestats(dist = LogNormal(), n = 40960)
```

```julia
@time plot_samplestats(dist = LogNormal(), n = 81920)
```

```julia
@time plot_samplestats(dist = LogNormal(), n = 163840)
```

```julia

```
