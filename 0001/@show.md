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
y = 2π
```

```julia
@show y
```

```julia
@which @show y
```

```julia
macro my_show(args...)
    A = [:(println($(string(x)), " = ", repr($(esc(x))))) for x in args]
    quote
        $(A...)
    end
end
```

```julia
@my_show y
```

```julia
x = "foo"
```

```julia
@show x y;
```

```julia
@my_show x y
```

```julia
@show sinpi(1/6);
```

```julia
@my_show sinpi(1/6)
```

```julia

```
