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
function name(dist::MixtureModel)
    c = components(dist)
    p = probs(dist)
    s = string(p[1]) * name(c[1])
    for i in 2:length(p)
        s = s * "+" * string(p[i]) * name(c[i])
    end
    s
end

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
function mcsimZW(;
        dist = Beta(0.2, 0.3),
        n = 20,
        L = 10^6,
    )
    μ, σ = mean(dist), std(dist)
    Z = Vector{Float64}(undef, L) # expected to follow Normal(1,0)
    W = similar(Z)
    tmp = [Vector{Float64}(undef, n) for i in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, tmp[Threads.threadid()])
        X̄ = mean(X) # sample mean
        S² = var(X) # sample unbiased variance
        Z[i] = √n * (X̄ - μ) / σ
        W[i] = √n * (S²/σ^2 - 1)
    end
    (; dist, n, Z, W)
end

function check_varcov(;
        dist = Beta(0.2, 0.3),
        n = 5,
        L = 10^8,
    )
    (; Z, W) = mcsimZW(; dist, n, L)
    Z² = @. Z^2
    κ₃ = myskewness(dist)
    κ₄ = mykurtosis(dist)
    
    @show dist
    @show κ₃
    @show κ₄
    meanZ,  meanZ_exact  = mean(Z),    0.0
    meanW,  meanW_exact  = mean(W),    0.0
    varZ,   varZ_exact   = var(Z),     1.0
    covZW,  covZW_exact  = cov(Z, W),  κ₃
    varW,   varW_exact   = var(W),     κ₄ + 2n/(n-1)
    covZ²W, covZ²W_exact = cov(Z², W), κ₄/√n
    varZ²,  varZ²_exact  = var(Z²),    κ₄/n + 2
    result = [
        :meanZ   meanZ   meanZ_exact   meanZ  - meanZ_exact
        :meanW   meanW   meanW_exact   meanW  - meanW_exact
        :varZ    varZ    varZ_exact    varZ   - varZ_exact
        :covZW   covZW   covZW_exact   covZW  - covZW_exact
        :varW    varW    varW_exact    varW   - varW_exact
        :covZ²W  covZ²W  covZ²W_exact  covZ²W - covZ²W_exact
        :varZ²   varZ²   varZ²_exact   varZ²  - varZ²_exact
    ]
end
```

```julia
check_varcov(dist = Normal(1, 2), n = 5)
```

```julia
check_varcov(dist = Beta(0.2, 0.3), n = 5)
```

```julia
check_varcov(dist = TDist(4.5), n = 5)
```

```julia
check_varcov(dist = TDist(4.6), n = 5)
```

```julia
check_varcov(dist = TDist(5.0), n = 5)
```

```julia
check_varcov(dist = Laplace(), n = 5)
```

```julia
check_varcov(dist = Exponential(), n = 5)
```

```julia
check_varcov(dist = Gamma(2, 1), n = 5)
```

```julia
check_varcov(dist = LogNormal(), n = 5)
```

```julia

```
