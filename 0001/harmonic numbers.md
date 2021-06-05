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

H_BigFloat(n) = harmonic(big(n))

function f(x, H=H_BigFloat)
    maxint = round(BigInt, exp(big(precision(BigFloat))))
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
result1 = [f(n) for n in 2:100]
```

```julia
using SpecialFunctions
H_BigFloat2(n) = digamma(big(n+1)) + MathConstants.γ
result2 = [f(n, H_BigFloat2) for n in 2:100]
```

```julia
result1 == result2
```

```julia
H_Float64(n) = harmonic(Float64(n))
result3 = [f(n, H_Float64) for n in 2:34]
```

```julia
A = getproperty.(result1[eachindex(result3)], :n)
B = getproperty.(result3, :n)
@show A[1:end-1] == B[1:end-1]
[(x = result1[k].x, n_result1 = A[k], n_result_2 = B[k]) for k in eachindex(A)]
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
