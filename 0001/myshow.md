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
a, b, c = 1, 2.0, "3"
@show a b c;
```

```julia
x = (a = 1, b = 2.0, c = "3")
```

```julia
function myshow(x::NamedTuple{names}) where names
    for name in names
        println(string(name), " = ", repr(getproperty(x, name)))
    end
end
```

```julia
myshow(x)
```

```julia

```
