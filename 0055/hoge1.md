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
    display_name: Julia
    language: julia
    name: julia
---

```julia
ENV["LINES"] = 200
ENV["COLUMNS"] = 200
using Printf
using Distributions
using LaTeXStrings
using StatsPlots
default(fmt=:png)
using SymPy
```

```julia
rd(x) = round(x; sigdigits=2)
n, ps, ks = 20, 0:0.1:1, 0:n
mat_doubled_cdf = @. rd(2cdf(Binomial(n, ps'), ks))
Any[
    "k\\p" permutedims(["$p" for p in ps])
    ks mat_doubled_cdf
]
```

```julia
rd(x) = round(x; sigdigits=4)
n, ps, ks = 20, 0:0.1:1, 0:n
mat_doubled_ccdf = @. rd(ccdf(Binomial(n, ps'), ks-1))
Any[
    "k\\p" permutedims(["$p" for p in ps])
    ks mat_doubled_ccdf
]
```

```julia
x ⪅ y = x < y || x ≈ y
x ⪆ y = x > y || x ≈ y

"""
    central法による二項分布の両側P値
二項分布binにおける「k以下になる確率の2倍」と「k以上になる確率の2倍」と「1」の中での最小の値

このP値からClopper-Pearsonの信頼区間が得られる。
"""
pvalue_central(bin, k) = min(1, 2cdf(bin, k), 2ccdf(bin, k-1))

"""
    minimum likelihood法による二項分布の両側P値
このP値からSterneの信頼区間が得られる。
"""
pvalue_minlike(bin, k) = sum(pdf(bin, i) for i in 0:n if pdf(bin, i) ⪅ pdf(bin, k))

"""
    score法による二項分布の両側P値
このP値は二項分布の正規分布近似で定義された両側P値と同じ。
"""
function pvalue_score(bin, k)
    (; n, p) = bin
    phat = k/n
    se = std(bin)
    z = (phat - p) / se
    2ccdf(Normal(), abs(z))
end
```

```julia
n, p, k = 20, 0.45, 14
pvalue_central(Binomial(n, p), k)
```

```julia
(; n, p) = bin
n, p
```

```julia
?pvalue_binomial_central
```

```julia
n, p = 20, 0.6
bin = Binomial(n, p)
bar(i -> pdf(bin, i), 8:n; label="")
bar!(i -> pdf(bin, i), 0:7; label="")
title!("Binomial($n, $p)")
```

```julia
2cdf(bin, 7)
```

```julia
2(0.03 + 0.13 + 0.49 + 1.46)/100
```

```julia
sum(pdf(bin, i) for i in 0:n if pdf(bin, i) ⪅ pdf(bin, 7))
```

```julia
normal = Normal(mean(bin), std(bin))
2cdf(normal, 7)
```

```julia
2cdf(normal, 7+0.5)
```

```julia
@show pvalue_central(bin, 7)
@show pvalue_minlike(bin, 7)
@show pvalue_score(bin, 7)
;
```

ほげ[^1]。


[^1] もげ
