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
code1 = raw"""
s = 0
for i in 1:10
    global s += i^3
end
println(s)
"""

include_string(Main, code1)
```

```julia
code2 = raw"""
s = 0
for i in 1:10
    s += i^3
end
println(s)
"""

include_string(Main, code2)
```

```julia
module A
s = 0
for i in 1:10
    global s += i^3
end
println(s)
end;
```

```julia
module B
s = 0
for i in 1:10
    s += i^3
end
println(s)
end;
```

```julia
s = 0

function f(n)
    for i in 1:n
        global s += i^3
    end
    println(s)
end

f(10)
```

```julia
function f(n)
    s = 0
    for i in 1:n
        s += i^3
    end
    println(s)
end

f(10)
```

```julia
s = 0

function f(n)
    for i in 1:n
        s += i^3
    end
    println(s)
end

f(10)
```

```julia
s = 0
let
    for i in 1:10
        global s += i^3
    end
    println(s)
end
```

```julia
let
    s = 0
    for i in 1:10
        s += i^3
    end
    println(s)
end
```

```julia
s = 0

let
    for i in 1:10
        s += i^3
    end
    println(s)
end
```

```julia

```
