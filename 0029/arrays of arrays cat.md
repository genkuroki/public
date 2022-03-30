---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Julia 1.8.0-beta1
    language: julia
    name: julia-1.8
---

```julia
@show VERSION

using BenchmarkTools

V = [rand(2, 3) for _ in 1:1000]

C1 = @btime cat($V...; dims=3)

C2 = @btime foldl((x, y) -> cat(x, y; dims=3), $V)

function aacat!(C, A::AbstractVector)
    for i in keys(A)
        C[axes(A[i])..., i] .= A[i]
    end
    C
end
function aacat!(C, A)
    for i in keys(A)
        C[axes(A[i])..., i.I...] .= A[i]
    end
    C
end
function aacat(A; A1 = A[begin])
    C = similar(A1, axes(A1)..., axes(A)...)
    aacat!(C, A)
end
C3 = @btime aacat($V)

C4 = similar(C3)
@btime aacat!($C4, $V)
aacat!(C4, V)

@show axes(C4)
@show C1 == C2 == C3 == C4;
```

```julia
using OffsetArrays
w = [OffsetArray(rand(2, 3), 0:1, -1:1) for _ in Iterators.product(1:3, 1:3)]
W = OffsetArray(w, 0:2, -2:0)
```

```julia
E = similar(W[begin], axes(W[begin])..., axes(W)...)
@btime aacat!($E, $W)
```

```julia
@btime aacat($W)
```

```julia

```
