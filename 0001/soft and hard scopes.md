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
# soft scope in REPL, Jupyter, etc.

n = 10
s = 0
for k in 1:n
    s += k^3
end
@show s;
```

```julia
# hard scope

module O1

n = 10
s = 0
for k in 1:n
    s += k^3
end
@show s

end;
```

```julia
# hard scope

module O2

n = 10
s = 0
for k in 1:n
    global s += k^3
end
@show s

end;
```

```julia
# hard scope

module O3

function f(n)
    s = 0
    for k in 1:n
        s += k^3
    end
    s
end

end

@show O3.f(10);
```

```julia

```
