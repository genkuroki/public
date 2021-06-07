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
A = [1 2 3 4]
```

```julia
B = [1, 2, 3, 4]
```

```julia
A*B
```

```julia
A = [1, 2, 3, 4]'
```

```julia
A*B
```

```julia
C = [
    1 5
    2 6
    3 7
    4 8
]
```

```julia
A*C
```

```julia
v = [
    1
    2
]
```

```julia
A*C*v
```

```julia
x = 1:4
y = 1:3
```

```julia
X = [x for y in y, x in x]
```

```julia
Y = [y for y in y, x in x]
```

```julia
zero(X)
```

```julia
module O

struct Foo{T}
    d::Dict{T}
end

function f(foo::Foo)
    d = foo.d
    for k in keys(d)
        println(k, ": ", d[k])
    end
end

end
```

```julia
foo = O.Foo(Dict(:a => 1, :b => 2, :c => 3))
```

```julia
O.f(foo)
```

```julia
@code_warntype O.f(foo)
```

```julia
module Q

struct Bar{K, V}
    d::Dict{K, V}
end

function g(bar::Bar)
    d = bar.d
    for k in keys(d)
        println(k, ": ", d[k])
    end
end

end
```

```julia
bar = Q.Bar(Dict(:a => 1, :b => 2, :c => 3))
```

```julia
Q.g(bar)
```

```julia
@code_warntype Q.g(bar)
```

```julia
bar.d[:d] = 4
```

```julia
bar
```

```julia
baz = (d = bar.d,)
```

```julia
baz.d[:e] = 5
baz
```

```julia
bar
```

```julia
struct Qoo{T} a::T end
qoo = Qoo(Dict{Symbol, Int64}())
```

```julia
qoo.a[:a] = 1
qoo
```

```julia
qoo.a[:b] = 2
qoo
```

```julia
qux = (a = Int64[],)
```

```julia
push!(qux.a, 1, 2, 3)
qux
```

```julia
qux.a[2] = 999
qux
```

```julia
f(x) = if x < 0; "negative"; elseif x == 0; "zero"; else; "positive"; end
```

```julia
f(-1), f(0), f(1)
```

```julia
g(x) = if x < 0 "negative" elseif x == 0 "zero" else "positive" end
```

```julia
g(-1), f(0), f(1)
```

```julia
:(if a A elseif b B elseif c C else D end) |> Base.remove_linenums!
```

```julia
hoe = (a = 1, b = 2.0, c = "three")
```

```julia
typeof(hoe)
```

```julia

```
