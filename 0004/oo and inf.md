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

<!-- #region tags=[] -->
## Python
<!-- #endregion -->

```julia tags=[]
from sympy import *
```

```julia
a = -oo
a
```

```julia
b = -float('inf')
b
```

```julia
type(a)
```

```julia
type(b)
```

## Julia

```julia
using SymPy
a = -oo
```

```julia
using PyCall
b = py"""-float('inf')"""
```

```julia
typeof(a)
```

```julia
typeof(b)
```

```julia

```
