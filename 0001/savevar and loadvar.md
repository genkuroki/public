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
savevar(fn, x) = write(fn, string(x))
loadvar(fn) = read(fn, String) |> Meta.parse |> eval
```

```julia
A = randn(ComplexF64, 4, 3, 2)
```

```julia
savevar("tmp/A.txt", A)
```

```julia
read("tmp/A.txt", String)
```

```julia
A_load = loadvar("tmp/A.txt")
```

```julia
A_load == A
```

```julia
B = ["Foo", "Bar", "Baz"]
```

```julia
savevar("tmp/B.txt", B)
```

```julia
read("tmp/B.txt", String)
```

```julia
B_load = loadvar("tmp/B.txt")
```

```julia
B_load == B
```

```julia
D = Dict(:A => A, :B => B)
```

```julia
savevar("tmp/D.txt", D)
```

```julia
read("tmp/D.txt", String)
```

```julia
D_load = loadvar("tmp/D.txt")
```

```julia
D_load == D
```

```julia
module O
struct Foo{A, B}
    a::A
    b::B
end
end
```

```julia
foo = O.Foo(A, B)
```

```julia
savevar("tmp/foo.txt", foo)
```

```julia
read("tmp/foo.txt", String)
```

```julia
foo_load = loadvar("tmp/foo.txt")
```

```julia
(foo_load.a, foo_load.b) == (foo.a, foo.b)
```

```julia

```
