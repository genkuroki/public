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
using Random
```

```julia
function simpi(N, rng = Random.default_rng())
    c = 0
    for _ in 1:N
        c += isone(gcd(rand(rng, Int), rand(rng, Int)))
    end
    √(6N/c)
end

simpi(10)
Random.seed!(4549373)

@time simpi(10^8)
```

```julia
Random.seed!(4549373)

@time let
    N = 10^8
    rng = Random.default_rng()
    s = 0
    for _ in 1:N
        s += isone(gcd(rand(rng, Int), rand(rng, Int)))
    end
    √(6N/s)
end
```

```julia
Random.seed!(4549373)

@time begin
    N = 10^8
    rng = Random.default_rng()
    s = 0
    for _ in 1:N
        s += isone(gcd(rand(rng, Int), rand(rng, Int)))
    end
    √(6N/s)
end
```

```julia
@show Threads.nthreads()

using Random

# https://github.com/genkuroki/MyUtils.jl
using MyUtils: @my_threads

function simpi_threads(N)
    a = Threads.Atomic{Int}(0)
    @my_threads begin
        rng = Random.default_rng()
        c = 0
    end for _ in 1:N
        c += isone(gcd(rand(rng, Int), rand(rng, Int)))
    end begin
        Threads.atomic_add!(a, c)
    end
    √(6N/a[])
end

simpi_threads(10)

@time simpi_threads(10^8)
```

```julia
@show Threads.nthreads()

using Random

# https://github.com/genkuroki/MyUtils.jl
using MyUtils: @my_threads

@time let
    N = 10^8
    a = Threads.Atomic{Int}(0)
    @my_threads begin
        rng = Random.default_rng()
        c = 0
    end for _ in 1:N
        c += isone(gcd(rand(rng, Int), rand(rng, Int)))
    end begin
        Threads.atomic_add!(a, c)
    end
    √(6N/a[])
end
```

```julia
@show Threads.nthreads()

using Random

# https://github.com/genkuroki/MyUtils.jl
using MyUtils: @my_threads

@time begin
    N = 10^8
    a = Threads.Atomic{Int}(0)
    @my_threads begin
        rng = Random.default_rng()
        c = 0
    end for _ in 1:N
        c += isone(gcd(rand(rng, Int), rand(rng, Int)))
    end begin
        Threads.atomic_add!(a, c)
    end
    √(6N/a[])
end
```

```julia

```
