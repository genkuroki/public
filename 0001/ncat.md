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
VERSION
```

```julia
using Random
using MetaUtils
```

```julia
A = [randstring(8) for i in 1:4, j in 1:3]
```

```julia
A_string = string(A)
print(A_string)
```

```julia
A_parse = Meta.parse(A_string)
```

```julia
show_expr(A_parse)
```

```julia
A_eval = eval(A_parse)
```

```julia
B = [randstring(8) for i in 1:4, j in 1:3, k in 1:2]
```

```julia
B_string = string(B)
print(B_string)
```

```julia
B_parse = Meta.parse(B_string)
```

```julia
using MetaUtils
show_expr(B_parse)
```

```julia
eval(B_parse)
```

```julia
VERSION
```

```julia
a = [
    111 112
    121 122
    ;;;
    211 212
    221 222
]
```

```julia
b = [
    "a" "b"
    "c" "d"
    ;;;
    "e" "f"
    "g" "h"
]
```

```julia
Base.ndims(x::AbstractString) = 0

b = [
    "a" "b"
    "c" "d"
    ;;;
    "e" "f"
    "g" "h"
]
```

```julia
c = [
    1       2.0
    "3"     '4'
    ;;;
    Int64   nothing
    missing sin
]
```

```julia

```
