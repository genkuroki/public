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
methods(&)
```

```julia
&(true, true, true, true)
```

```julia
Base.:&(true, true, true, true)
```

```julia
var"&"(true, true, true, true)
```

```julia
.&(true, true, true, true)
```

```julia
.&(true, true, false, true)
```

```julia
@code_llvm debuginfo=:none Base.:&(true, true, true, true)
```

```julia
f(x...) = .&(x...)
@code_llvm debuginfo=:none f(true, true, true, true)
```

```julia

```
