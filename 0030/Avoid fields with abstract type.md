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
module O

"""
Implement `elzero(x)`, `Base.:+(x, y)`, `Base.length(x)`, and `Base.getindex(x, i)` for a subtype `T` of `AbstractFoo`.  Then you can use `naivesum(x)` method for `x::T`.
"""
abstract type AbstractFoo end
function naivesum(x::AbstractFoo)
    s = elzero(x)
    for i in 1:length(x)
        s += x[i]
    end
    s
end

struct Foo <: AbstractFoo
    a
end
elzero(x::Foo) = Foo(zero(eltype(x.a)))
Base.:+(x::Foo, y::Foo) = Foo(x.a + y.a)
Base.length(x::Foo)::Int = length(x.a)
Base.getindex(x::Foo, i::Integer) = Foo(x.a[i])

struct Bar{T} <: AbstractFoo
    a::T
end
elzero(x::Bar) = Bar(zero(eltype(x.a)))
Base.:+(x::Bar, y::Bar) = Bar(x.a + y.a)
Base.length(x::Bar) = length(x.a)
Base.getindex(x::Bar, i::Integer) = Bar(x.a[i])

end

@doc O.AbstractFoo
```

```julia
a = randn(10^6)
foo = O.Foo(a)
bar = O.Bar(a);
```

```julia
@code_warntype O.naivesum(foo)
```

```julia tags=[]
@code_warntype O.naivesum(bar)
```

```julia
@time O.naivesum(foo)
@time O.naivesum(foo)
@time O.naivesum(foo)
```

```julia
@time O.naivesum(bar)
@time O.naivesum(bar)
@time O.naivesum(bar)
```

```julia

```
