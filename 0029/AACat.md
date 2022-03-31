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
module O

struct AACat{T, N, U} <: AbstractArray{T, N}
    A::U
end
AACat(A) = AACat{eltype(eltype(A)), ndims(eltype(A)) + ndims(A), typeof(A)}(A)

function Base.size(C::AACat)
    A = C.A
    (size(A[begin])..., size(A)...)
end

function Base.axes(C::AACat)
    A = C.A
    (axes(A[begin])..., axes(A)...)
end

function Base.getindex(C::AACat, I::Integer...)
    A = C.A
    n = length(I) - ndims(A)
    J = I[1:n]
    K = I[n+1:end]
    getindex(A[K...], J...)
end

function Base.setindex!(C::AACat, v, I::Integer...)
    A = C.A
    n = length(I) - ndims(A)
    J = I[1:n]
    K = I[n+1:end]
    setindex!(A[K...], v, J...)
end

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

end
```

```julia
A = [rand(0:9, 2, 3) for _ in Iterators.product(1:3, 1:2)]
```

```julia
AC = O.AACat(A)
```

```julia
using OffsetArrays

b = [OffsetArray(rand(0:9, 2, 3), 0:1, 0:2) for _ in Iterators.product(1:3, 1:2)]
B = OffsetArray(b, 0:2, 0:1)
```

```julia
BC = O.AACat(B)
```

```julia
BC[1, 1, :, :]
```

```julia
BC[:, :, 1, 1]
```

```julia
BC[1, 1, 1, 1] = 99
BC[:, :, 1, 1]
```

```julia
using SplitApplyCombine

c = [OffsetArray(rand(0:9, 2, 3), 0:1, 0:2) for _ in Iterators.product(1:3, 1:2)]
C = OffsetArray(c, 0:2, 0:1)

CC = combinedims(C)
```

```julia
CC[1, 1, :, :]
```

```julia
CC[:, :, 1, 1]
```

```julia
CC[1, 1, 1, 1] = 99
CC[:, :, 1, 1]
```

```julia
d = [OffsetArray(rand(0:9, 2, 3), 0:1, 0:2) for _ in Iterators.product(1:3, 1:2)]
D = OffsetArray(c, 0:2, 0:1)

DC = combinedimsview(D)
```

```julia
DC[1, 1, :, :]
```

```julia
DC[:, :, 1, 1]
```

```julia
DC[1, 1, 1, 1] = 99
```

```julia
using RecursiveArrayTools

V = [rand(0:9, 2, 3) for _ in 1:4]
VC = VectorOfArray(V)
```

```julia
[
    VC[1,1,1] VC[1,2,1] VC[1,3,1]
    VC[2,1,1] VC[2,2,1] VC[2,3,1]
]
```

```julia
using BenchmarkTools

V = [rand(2, 3) for _ in 1:1000]

println("---------- O.aacat")
C3 = @btime O.aacat($V)
s3 = @btime sum($C3)
@show s3

println("---------- O.AACat")
C5 = @btime O.AACat($V)
s5 = @btime sum($C5)
@show s5

println("---------- SplitApplyCombine.combinedims")
C6 = @btime combinedims($V)
s6 = @btime sum($C6)
@show s6

println("---------- SplitApplyCombine.combinedimsview")
C7 = @btime combinedimsview($V)
s7 = @btime sum($C7)
@show s7

println("---------- RecursiveArrayTools.VectorOfArray")
C8 = @btime VectorOfArray($V)
s8 = @btime sum($C8)
@show s8

println("----------")
@show s3 ≈ s5 ≈ s6 ≈ s7 ≈ s8;
```

```julia
using BenchmarkTools

V = [rand(2, 3) for _ in 1:10^6]

println("---------- O.aacat")
C3 = @btime O.aacat($V)
s3 = @btime sum($C3)
@show s3

println("---------- O.AACat")
C5 = @btime O.AACat($V)
s5 = @btime sum($C5)
@show s5

println("---------- SplitApplyCombine.combinedims")
C6 = @btime combinedims($V)
s6 = @btime sum($C6)
@show s6

println("---------- SplitApplyCombine.combinedimsview")
C7 = @btime combinedimsview($V)
s7 = @btime sum($C7)
@show s7

println("---------- RecursiveArrayTools.VectorOfArray")
C8 = @btime VectorOfArray($V)
s8 = @btime sum($C8)
@show s8

println("----------")
@show s3 ≈ s5 ≈ s6 ≈ s7 ≈ s8;
```

```julia

```
