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
struct Foo end
```

```julia tags=[]
Base.getproperty(::Type{Foo}, f::Symbol) = println("Foo's property: ", f)
```

```julia
Foo.qoo
```

```julia
Union{}.qoo
```

```julia
struct Bar end
```

```julia
Base.getproperty(::Type{T}, f::Symbol) where T<:Bar = println("Bar's property: ", f)
```

```julia
Bar.qoo
```

```julia
Union{}.qoo
```

```julia
Union{} <: Bar
```

__type piracy!__

See https://github.com/JuliaLang/julia/issues/39534 and https://github.com/JuliaLang/julia/pull/39573

```julia

```
