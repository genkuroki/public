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
expr = :(if a A elseif b B elseif c C else D end) |> Base.remove_linenums!
```

```julia
dump(expr)
```

```julia
using MetaUtils
```

```julia
MetaUtils.print_tree(expr)
```

```julia
show_expr(expr)
```

```julia

```
