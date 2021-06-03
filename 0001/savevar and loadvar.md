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

const dir_savevar = "tmp"
fn_savevar(x::Symbol) = joinpath(dir_savevar, string(x) * ".txt")
macro savevar(x) :(savevar($(fn_savevar(x)), $(esc(x)))) end
macro loadvar(x) :(loadvar($(fn_savevar(x)))) end
```

```julia
A = randn(ComplexF64, 4, 3, 2)
```

```julia
@savevar A
```

```julia
read("tmp/A.txt", String)
```

```julia
A_load = @loadvar A
```

```julia
A_load == A
```

```julia
B = ["Foo", "Bar", "Baz"]
```

```julia
@savevar B
```

```julia
read("tmp/B.txt", String)
```

```julia
B_load = @loadvar B
```

```julia
B_load == B
```

```julia
D = Dict(:A => A, :B => B)
```

```julia
@savevar D
```

```julia
read("tmp/D.txt", String)
```

```julia
D_load = @loadvar D
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
@savevar foo
```

```julia
read("tmp/foo.txt", String)
```

```julia
foo_load = @loadvar foo
```

```julia
(foo_load.a, foo_load.b) == (foo.a, foo.b)
```

```julia

```
