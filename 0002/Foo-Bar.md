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
struct Foo{A} a::A end
foo = Foo(Int[])
push!(foo.a, 123)
foo
```

```julia
push!(foo.a, 456)
foo.a[1] = 999
foo
```

```julia
Base.@kwdef struct Bar{A, B, C}
     a::A = 1
     b::B = 2.0
     c::C = "three"
end
bar = Bar()
```

```julia
(; c, a) = bar
@show c a;
```

```julia

```
