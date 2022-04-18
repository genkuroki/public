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

* Sterne, Theodore E. Some Remarks on Confidence or Fiducial Limits. Biometrika
Vol. 41, No. 1/2 (Jun., 1954), pp. 275-278 (4 pages) https://www.jstor.org/stable/2333026

```julia
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10)
```

```julia
function pvalue_clopper_pearson(dist::DiscreteUnivariateDistribution, x)
    min(1, 2cdf(dist, x), 2ccdf(dist, x-1))
end
pvalue_clopper_pearson(n, k, p) = pvalue_clopper_pearson(Binomial(n, p), k)
```

```julia
x ⪅ y = x < y || x ≈ y

# Naive implementation is terribly slow.
function pvalue_stern_naive(dist::DiscreteUnivariateDistribution, x; xmax = 10^6)
    Px = pdf(dist, x)
    Px == 0 && return Px
    ymin, maxdist = minimum(dist), maximum(dist)
    ymax = maxdist == Inf ? xmax : maxdist
    sum(pdf(dist, y) for y in ymin:ymax if 0 < pdf(dist, y) ⪅ Px; init = 0.0)
end
pvalue_stern_naive(n, k, p) = pvalue_stern_naive(Binomial(n, p), k)

# Second implementation is very slow.
function pvalue_stern_old(dist::DiscreteUnivariateDistribution, x)
    Px = pdf(dist, x)
    Px == 0 && return Px
    distmin, distmax = extrema(dist)
    m = mode(dist)
    Px ≈ pdf(dist, m) && return one(Px)
    if x < m
        y = m + 1
        while !(pdf(dist, y) ⪅ Px)
            y += 1
        end
        cdf(dist, x) + ccdf(dist, y-1)
    else # k > m
        y = m - 1
        while !(pdf(dist, y) ⪅ Px)
            y -= 1
        end
        cdf(dist, y) + ccdf(dist, x-1)
    end
end
pvalue_stern_old(n, k, p) = pvalue_stern_old(Binomial(n, p), k)

### efficient implementation

_pdf_le(x, (dist, y)) =  pdf(dist, x) ⪅ y

function _search_boundary(f, x0, Δx, param)
    x = x0
    if f(x, param)
        while f(x - Δx, param) x -= Δx end
    else
        x += Δx
        while !f(x, param) x += Δx end
    end
    x
end

function pvalue_stern(dist::DiscreteUnivariateDistribution, x)
    Px = pdf(dist, x)
    Px == 0 && return Px
    m = mode(dist)
    Px ≈ pdf(dist, m) && return one(Px)
    if x < m
        y = _search_boundary(_pdf_le, 2m - x, 1, (dist, Px))
        cdf(dist, x) + ccdf(dist, y-1)
    else # x > m
        y = _search_boundary(_pdf_le, 2m - x, -1, (dist, Px))
        cdf(dist, y) + ccdf(dist, x-1)
    end
end
pvalue_stern(n, k, p) = pvalue_stern(Binomial(n, p), k)
```

```julia
n = 10
k = -1:11
p = 0.4
a = @time pvalue_stern_naive.(n, k, p)
b = @time pvalue_stern_old.(n, k, p)
c = @time pvalue_stern.(n, k, p)
d = @time pvalue_clopper_pearson.(n, k, p)
@show a ≈ b ≈ c
[a b c d]
```

```julia
n = 100000
k = 49500:50500
a = @time pvalue_stern_naive.(n, k, 0.5)
b = @time pvalue_stern_old.(n, k, 0.5)
c = @time pvalue_stern.(n, k, 0.5)
d = @time pvalue_clopper_pearson.(n, k, 0.5)
@show a ≈ b ≈ c ≈ d;
```

```julia
dist = Hypergeometric(9, 9, 9)
ran = -1:10
a = @time pvalue_stern_naive.(dist, ran)
b = @time pvalue_stern_old.(dist, ran)
c = @time pvalue_stern.(dist, ran)
d = @time pvalue_clopper_pearson.(dist, ran)
@show a ≈ b ≈ c ≈ d
[a b c d]
```

```julia
dist = Poisson(4)
ran = -1:10
a = @time pvalue_stern_naive.(dist, ran)
c = @time pvalue_stern.(dist, ran)
d = @time pvalue_clopper_pearson.(dist, ran)
[a c d]
```

```julia

```
