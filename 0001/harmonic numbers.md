---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Julia 1.7.0-DEV
    language: julia
    name: julia-1.7
---

```julia
using Polylogarithms
setprecision(200)

H(n) = harmonic(big(n))

function f(x)
    maxint = round(BigInt, exp(big(200)))
    a, b = one(BigInt), maxint
    @assert H(a) < x
    @assert H(b) ≥ x
    maxiters = 2floor(Int, log2(maxint))
    for _ in 1:maxiters
        a + 1 == b && break
        c = a + (b - a) ÷ 2
        if H(a) < x ≤ H(c)
            b = c
        else
            a = c
        end
    end
    (x = x, n = b, harmonic_num = H(b))
end

ENV["LINES"] = 256
[f(n) for n in 2:100]
```

```julia
using SpecialFunctions

harmonic_naive(n) = sum(inv, Base.OneTo(n))
harmonic_digamma(n) = digamma(n + 1) - digamma(1)
```

```julia
harmonic_naive(10)
```

```julia
harmonic_digamma(10)
```

```julia
[(n, harmonic_digamma(n) - harmonic_naive(n)) for n in 1:20]
```

```julia
[(2^n, harmonic_digamma(2^n) - harmonic_naive(2^n)) for n in 1:20]
```

```julia
using Base.MathConstants: γ
using Plots

harmonic_approx(x) = log(x) + γ + 1/(2x)

n = 1:50
x = range(extrema(n)...; length=1000)
plot(; legend=:bottomright)
scatter!(n, harmonic_naive.(n); label="harmonic numbers")
plot!(x, harmonic_approx.(x); label="log(x) + γ + 1/(2x)")
```

```julia
[(n, harmonic_approx(n) - harmonic_naive(n)) for n in 1:20]
```

```julia
harmonic_approx2(x) = log(x) + γ + 1/(2x) - 1/(12x^2)

[(n, harmonic_approx2(n) - harmonic_naive(n)) for n in 1:20]
```

```julia

```
