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

```
