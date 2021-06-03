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
module My

"""
    savevar(fn, x)

saves the value of `x` to the file `fn`, where `fn` is the filename string of the file.
"""
savevar(fn, x) = write(fn, string(x))

"""
    loadvar(fn)

loads the file `fn` (the filename string of the file) and `Meta.parse |> eval`.
"""
loadvar(fn) = read(fn, String) |> Meta.parse |> eval

"""
    dir_savevar[]

is the default directory to which `@savevar` saves the values of variables.
"""
const dir_savevar = Ref(".")

"""
    fn_savevar(x::Symbol)

is the filename string to which `@savevar` saves the value of a variable.
"""
fn_savevar(x::Symbol) = joinpath(dir_savevar[], string(x) * ".txt")

"""
    @savevar(args...)

saves the variables in args to the corresponding textfiles.

Example: `@savevar A B C` saves the variables `A`, `B`, `C` to textfiles. 
"""
macro savevar(args...)
    A = [:(savevar($(fn_savevar(x)), $(esc(x)))) for x in args]
    quote $(A...); nothing end
end

"""
    @loadvar(args...)

loads the values from the textfiles corresponding to `args`.
If `length(args)` is greater than 1, then it returns the tuple of the values.

Example: `a, b, c = @loadvar A B C` loads 
the values of `A`, `B`, `C` in textfiles to the variables `a`, `b`, `c`.
"""
macro loadvar(args...)
    if length(args) == 1
        x = args[1]
        :(loadvar($(fn_savevar(x))))
    else
        A = [:(loadvar($(fn_savevar(x)))) for x in args]
        :(($(A...),))
    end
end

end
```

```julia
using .My: dir_savevar, @savevar, @loadvar
dir_savevar[] = "tmp"
```

```julia
?dir_savevar
```

```julia
?@savevar
```

```julia
?@loadvar
```

```julia
using Random; Random.seed!(4649373)
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
x, y, z = randn(3)
x, y, z
```

```julia
@savevar x y z
```

```julia
X, Y, Z = @loadvar x y z
X, Y, Z
```

```julia
(X, Y, Z) == (x, y, z)
```

```julia

```
