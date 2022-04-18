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
using Printf
using Random
Random.seed!(4649373)
using StaticArrays
using StatsBase
using StatsPlots
default(fmt = :png, titlefontsize=10, size = (400, 250))
```

```julia
function chisqstat(x, n, p)
    sum((x[i] - n*p[i])^2/(n*p[i]) for i in eachindex(p))
end

pval(r, chi2) = ccdf(Chisq(r - 1), chi2)
pval(x, n, p) = pval(length(p), chisqstat(x, n, p))

function mcsim(n, p₀, p = p₀; L = 10^6)
    dist_true = Multinomial(n, p₀)
    C = Vector{Float64}(undef, L)
    P = similar(C)
    Z = Vector{Int}(undef, L)
    r = length(p₀)
    tmp = [Vector{Int}(undef, r) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist_true, tmp[Threads.threadid()])
        C[i] = chisqstat(X, n, p)
        P[i] = pval(r, C[i])
        Z[i] = count(==(0), X)
    end
    ecdfC, ecdfP = ecdf(C), ecdf(P)
    countZ = count.(.==(0:r-1), Ref(Z))
    ecdfC, ecdfP, countZ, C, P, Z
end

function plot_mcsim(n, p₀, p = p₀; L = 10^6)
    r = length(p₀)
    ecdfC, ecdfP, countZ, C, P, Z = mcsim(n, p₀, p; L)
    prob_numzero = countZ / L
    
    P1 = plot(; legend = :bottomright)
    plot!(x -> ecdfC(x), 0, r+4√(2r); label="ecdf(χ²-stats)")
    plot!(x -> cdf(Chisq(r-1), x); label="cdf(Chisq($(r-1)))", ls=:dash)

    P2 = plot(; legend = :topleft)
    plot!(x -> ecdfP(x), 0, 0.1; label="ecdf(p-values)")
    plot!(x -> cdf(Uniform(), x); label="cdf(Uniform())", ls=:dash)
    
    @show n
    @show p₀
    p != p₀ && @show p
    @show n*p₀
    p != p₀ && @show n*p
    @show prob_numzero
    plot(P1, P2; size=(600, 250))
end

eprob(v) = v / sum(v)

function evenprob(t, v = ones(6))
    s = sum(v)
    a, b, c = (v[1]+v[2], v[3]+v[4], v[5]+v[6]) ./ s
    SVector((1-t)*a, t*a, (1-t)*b, t*b, (1-t)*c, t*c)
end
```

```julia
plot_mcsim(10, evenprob(0.75))
```

```julia
plot_mcsim(20, evenprob(0.75))
```

```julia
plot_mcsim(30, evenprob(0.75))
```

```julia
plot_mcsim(100, evenprob(0.75))
```

```julia
plot_mcsim(24, evenprob(0.5))
```

```julia
plot_mcsim(10, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_mcsim(20, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_mcsim(30, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_mcsim(100, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_mcsim(100, eprob([2, 8]))
```

```julia
function plot_pvaluefunction!(x, n, p; kwargs...)
    plot!(f, 0, 1; label="", kwargs...)
end

function plot_pvaluefunctions(n, p₀; L=100)
    @show n
    @show p₀
    @show n*p₀
    dist_true = Multinomial(n, p₀)
    
    P1 = plot(; xtick=0:0.1:1, ytick=0:0.1:1)
    P2 = plot(; xtick=0:0.1:1, ytick=0:0.1:1)
    for _ in 1:L
        X = rand(dist_true)
        f(t) = pval(X, n, evenprob(t, p₀))
        plot!(P1, f, 0, 1; label="",  alpha=0.4, lw=0.5)
        g(t) = pval(SVector(X[1]+X[3]+X[5], X[2]+X[4]+X[6]), n, SVector(1-t, t))
        plot!(P2, g, 0, 1; label="",  alpha=0.4, lw=0.5)
    end
    title!(P1, "multinomial p-value functions (n = $n)")
    title!(P2, "binomial p-value functions (n = $n)")
    plot(P1, P2; size=(800, 250))
end
```

```julia
plot_pvaluefunctions(10, evenprob(0.75))
```

```julia
plot_pvaluefunctions(20, evenprob(0.75))
```

```julia
plot_pvaluefunctions(30, evenprob(0.75))
```

```julia
plot_pvaluefunctions(100, evenprob(0.75))
```

```julia
plot_pvaluefunctions(400, evenprob(0.75))
```

```julia
plot_pvaluefunctions(10, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_pvaluefunctions(20, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_pvaluefunctions(30, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_pvaluefunctions(100, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_pvaluefunctions(400, eprob([1, 2, 3, 5, 8, 13]))
```

```julia
plot_pvaluefunctions(10, evenprob(0.5))
```

```julia
plot_pvaluefunctions(20, evenprob(0.5))
```

```julia
plot_pvaluefunctions(30, evenprob(0.5))
```

```julia
plot_pvaluefunctions(100, evenprob(0.5))
```

```julia
plot_pvaluefunctions(400, evenprob(0.5))
```

```julia
function anim_pvaluefunctions(n, p₀; L=100)
    @show n
    @show p₀
    @show n*p₀
    dist_true = Multinomial(n, p₀)
    
    anim = @animate for _ in 1:L
        plot(; xtick=0:0.1:1, ytick=0:0.1:1)
        X = rand(dist_true)
        f(t) = pval(X, n, evenprob(t, p₀))
        plot!(f, 0, 1; label="mult")
        g(t) = pval(SVector(X[1]+X[3]+X[5], X[2]+X[4]+X[6]), n, SVector(1-t, t))
        plot!(g, 0, 1; label="bin")
    end
    
    gif(anim, "anim_pvaluefunctions.gif"; fps=2)
end
```

```julia
anim_pvaluefunctions(10, evenprob(0.5))
```

```julia
anim_pvaluefunctions(20, evenprob(0.5))
```

```julia
anim_pvaluefunctions(30, evenprob(0.5))
```

```julia
anim_pvaluefunctions(100, evenprob(0.5))
```

```julia
anim_pvaluefunctions(400, evenprob(0.5))
```

```julia
function f1()
    function g()
        x = 1
        return () -> x
    end
    x = 0
    h = g()
    return h()
end

f1()
```

```julia
function f2()
    function g()
        x = 1
        return () -> x
    end
    h = g()
    x = 0
    return h()
end

f2()
```

```julia
@code_typed f1()
```

```julia
@code_typed f2()
```

```julia

```
