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
x ⪅ y = x < y || x ≈ y

# Naive implementation is terrible slow.
function pvalue_stern_naive(dist::DiscreteUnivariateDistribution, x; xmax = 10^6)
    Px = pdf(dist, x)
    Px == 0 && return Px
    m = mode(dist)
    (x == m || Px ≈ pdf(dist, m)) && return 1.0
    ymin, maxdist = minimum(dist), maximum(dist)
    ymax = maxdist == Inf ? xmax : maxdist
    sum(pdf(dist, y) for y in ymin:ymax if 0 < pdf(dist, y) ⪅ Px; init = 0.0)
end

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
    (x == m || Px ≈ pdf(dist, m)) && return 1.0
    if x < m
        y = _search_boundary(_pdf_le, 2m - x, 1, (dist, Px))
        cdf(dist, x) + ccdf(dist, y-1)
    else # x > m
        y = _search_boundary(_pdf_le, 2m - x, -1, (dist, Px))
        cdf(dist, y) + ccdf(dist, x-1)
    end
end
```

```julia
dist = Binomial(10, 0.3)
a = pvalue_stern.(dist, -1:11)
b = pvalue_stern_naive.(dist, -1:11)
@show a ≈ b
a
```

```julia
dist = Binomial(100000, 0.5)
ran = 49500:50500
@time a = pvalue_stern.(dist, ran)
@time b = pvalue_stern_naive.(dist, ran)
@show a ≈ b
a
```

```julia
dist = Hypergeometric(10, 10, 10)
a = pvalue_stern.(dist, -1:11)
b = pvalue_stern_naive.(dist, -1:11)
@show a ≈ b
a
```

```julia
dist = Poisson(4)
a = pvalue_stern.(dist, -1:10)
b = pvalue_stern_naive.(dist, -1:10)
@show a ≈ b
a
```

```julia

```
