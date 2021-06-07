---
jupyter:
  jupytext:
    formats: ipynb,md,jl:hydrogen
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
module Q

const Foo = NamedTuple{(:a,)}
const Bar = NamedTuple{(:a, :b)}
const Baz = NamedTuple{(:a, :b, :c)}

f(x::Foo) = println("Foo", x)
f(x::Bar) = println("Bar", x)
f(x::Baz) = println("Baz", x)

end

foo = (a = 1,)
bar = (a = 1, b = 2.0)
baz = (a = 1, b = 2.0, c = "three");

Q.f(foo)
Q.f(bar)
Q.f(baz);
```

```julia
baz = (a = 999, baz.b, baz.c)
baz
```

```julia
qoo = (a = Ref(1), b = Ref(2.0), c = Ref("three"))
```

```julia
qoo.a[] = 1000
qoo.b[] = 2000.0
qoo.c[] = "three thousand"
qoo
```

```julia
qux = (a = fill(1), b = fill(2.0), c = fill("three"))
```

```julia
qux.a[] = 1000
qux.b[] = 2000.0
qux.c[] = "three thousand"
qux
```

```julia

```
