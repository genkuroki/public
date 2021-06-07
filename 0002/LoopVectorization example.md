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
Threads.nthreads()
```

```julia
using LoopVectorization
using BenchmarkTools

function L!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    for j in J
        for i in I
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function M!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    for i in I
        for j in J
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function L_inbounds!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    @inbounds for j in J
        for i in I
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function M_inbounds!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    @inbounds for i in I
        for j in J
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function L_turbo!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    @turbo for j in J
        for i in I
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function M_turbo!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    @turbo for i in I
        for j in J
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function L_turbo2!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    for j in J
        @turbo for i in I
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function M_turbo2!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    for i in I
        @turbo for j in J
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function L_tturbo!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    @tturbo for j in J
        for i in I
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function M_tturbo!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    @tturbo for i in I
        for j in J
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function L_tturbo2!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    for j in J
        @tturbo for i in I
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end

function M_tturbo2!(foo, bar)
    I, J = axes(foo, 1)[begin+1:end-1], axes(foo, 2)[begin+1:end-1]
    for i in I
        @tturbo for j in J
            foo[i, j] = bar[i, j-1] + bar[i, j+1] + bar[i-1, j] + bar[i+1, j] - 4bar[i, j]
        end
    end
end
```

```julia
n = 100
x = y = range(-1, 1; length=n)
f(x, y) = exp(-2(x^2+y^2))
bar = f.(x', y)

foo = zero(bar)
foo_inbounds = zero(bar)
foo_turbo = zero(bar)
foo_turbo2 = zero(bar)
foo_tturbo = zero(bar)
foo_tturbo2 = zero(bar)
Foo = zero(bar)
Foo_inbounds = zero(bar)
Foo_turbo = zero(bar)
Foo_turbo2 = zero(bar)
Foo_tturbo = zero(bar)
Foo_tturbo2 = zero(bar)

L!(foo, bar)
L_inbounds!(foo_inbounds, bar)
L_turbo!(foo_turbo, bar)
L_turbo2!(foo_turbo2, bar)
L_tturbo!(foo_tturbo, bar)
L_tturbo2!(foo_tturbo2, bar)
M!(Foo, bar)
M_inbounds!(Foo_inbounds, bar)
M_turbo!(Foo_turbo, bar)
M_turbo2!(Foo_turbo2, bar)
M_tturbo!(Foo_tturbo, bar)
M_tturbo2!(Foo_tturbo2, bar)

foo == foo_inbounds == foo_turbo == foo_turbo2 == foo_tturbo == foo_tturbo2 ==
Foo == Foo_inbounds == Foo_turbo == Foo_turbo2 == Foo_tturbo == Foo_tturbo2
```

```julia
@btime L!($foo, $bar)
```

```julia
@btime M!($foo, $bar)
```

```julia
@btime L_inbounds!($foo, $bar)
```

```julia
@btime M_inbounds!($foo, $bar)
```

```julia
@btime L_turbo!($foo, $bar)
```

```julia
@btime M_turbo!($foo, $bar)
```

```julia
@btime L_turbo2!($foo, $bar)
```

```julia
@btime M_turbo2!($foo, $bar)
```

```julia
@btime L_tturbo!($foo, $bar)
```

```julia
@btime M_tturbo!($foo, $bar)
```

```julia
@btime L_tturbo2!($foo, $bar)
```

```julia
@btime M_tturbo2!($foo, $bar)
```

```julia

```
