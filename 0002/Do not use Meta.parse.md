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

https://twitter.com/wasisama/status/1401804568468484099

```julia
ENV["LINES"] = 20
ENV["COLUMNS"] = 100

using MyUtils: showimg
showimg("image/png", "E3Q1019VoAEZfmL.png")
```

Don't use Meta.parse!

See [How to warn new users away from metaprogramming](https://discourse.julialang.org/t/how-to-warn-new-users-away-from-metaprogramming/35022).

```julia
function equalrows(a, r, names = propertynames(a))
    a[.&((getproperty(a, n) .== getproperty(r, n) for n in names)...), :]
end
```

```julia
using DataFrames

names = [:index, :aaa, :bbb, :ccc, :ddd, :eee, :wfawefa, :okfkpe, :eorjgr]

data_a = Vector[[collect(1:10)]; [rand(10) for _ in 2:length(names)]]
data_r = deepcopy(data_a)
data_a[2][1:2:10] .= -1.0
data_a[7][1:3:10] .= -2.0

a = DataFrame(data_a, names)
r = DataFrame(data_r, names);
```

```julia
a
```

```julia
r
```

```julia
@show names
equalrows(a, r)
```

```julia
@show names_aaa_to_eee = [:aaa, :bbb, :ccc, :ddd, :eee]
equalrows(a, r, names_aaa_to_eee)
```

```julia
@show names_last_3 = [:wfawefa, :okfkpe, :eorjgr]
equalrows(a, r, names_last_3)
```

```julia
function equalrows_revised(a, r, names = propertynames(a))
    a[.&((isequal.(getproperty(a, n), getproperty(r, n)) for n in names)...), :]
end
```

```julia
@show names_last_3 = [:wfawefa, :okfkpe, :eorjgr]
equalrows_revised(a, r, names_last_3)
```

```julia
A, R = copy(a), copy(r)
A.okfkpe[2] = R.okfkpe[2] = NaN;
```

```julia
A
```

```julia
R
```

```julia
equalrows(A, R, names_last_3)
```

```julia
equalrows_revised(A, R, names_last_3)
```

```julia
B, S = copy(a), copy(r)
tmp = Vector{Union{Missing, eltype(B.okfkpe)}}(B.okfkpe)
B.okfkpe = S.okfkpe = tmp
B.okfkpe[2] = S.okfkpe[2] = missing;
```

```julia
B
```

```julia
S
```

```julia
equalrows(B, S, names_last_3)
```

```julia
equalrows_revised(B, S, names_last_3)
```

```julia

```
