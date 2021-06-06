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

__References__

* https://github.com/genkuroki/public/blob/main/0002/julia%20vs.%20gcc%20-%202021-06-06%20harmonic%20number.ipynb
* https://twitter.com/genkuroki/status/1401330514175291396
* https://github.com/genkuroki/public/blob/main/0001/harmonic%20numbers.ipynb
* https://twitter.com/genkuroki/status/1400995381933051904

```julia
function add_kbn(s, c, a)
    t = s + a
    d = c + ifelse(abs(s) ≥ abs(a), (s-t) + a, (a-t) + s)
    t, d
end

function sum_kbn(f, iter; T = Float64, s = zero(T), c = zero(T))
    for n in iter
        s, c = add_kbn(s, c, f(T(n)))
    end
    s + c
end

sum_kbn(iter; T = Float64, s = zero(T), c = zero(T)) =
    sum_kbn(identity, iter; T, s, c)

function sumupto_kbn(f, iter, x; T = Float64, s = zero(T), c = zero(T))
    for n in iter
        s, c = add_kbn(s, c, f(T(n)))
        s + c ≥ x && return s + c, n
    end
    s + c, typemax(eltype(iter))
end

@time sum_kbn(inv, 1:10^9)
```

```julia
using SpecialFunctions, Printf
H(n) = digamma(n+1) + MathConstants.γ
@printf "%.15f" H(big(10^9))
```

```julia
@time s = sum_kbn(inv, 1:7*10^8)
@time h, n = sumupto_kbn(inv, 7*10^8+1:10^9, 21.0; s)
```

```julia
@printf "%.15f" H(big(n))
```

```julia
@show Threads.nthreads()
using Distributed

function sum_kbn_threads(f, ran; T = Float64, nth = Threads.nthreads())
    splitran = Distributed.splitrange(ran[begin], ran[end], nth)
    S, C = zeros(T, nth), zeros(T, nth)
    Threads.@threads for i in 1:nth
        S[i] = sum_kbn(f, splitran[i]; T)
    end
    S .+= C
    S
end
```

```julia
@time S = sum_kbn_threads(inv, 1:10^9)
sum_kbn(S)
```

```julia
m, x = 7*10^8, 21.0
@time S = sum_kbn_threads(inv, 1:m)
@time s = sum_kbn(S)
@time h, n = sumupto_kbn(inv, m+1:round(Int, 1.1*m), x; s)
```

```julia
@printf "%.15f" H(big(n))
```

```julia
@time S = sum_kbn_threads(inv, 1:6*10^9)
sum_kbn(S)
```

```julia
m, x = 6*10^12, 30.0
@time S = sum_kbn_threads(inv, 1:m)
@time s = sum_kbn(S)
@time h, n = sumupto_kbn(inv, m+1:round(Int, 1.1*m), x; s)
```

```julia
@printf "%.15f" H(big(n))
```

```julia

```
