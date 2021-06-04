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
num_list = [1, 2, 3, 4, 5]
```

`[1, 2, 3, 4, 5]` はリストではなく、配列(もしくはベクトル)である。

```julia
num_arr = [1, 2, 3, 4, 5]
```

```julia
function f(num_arr)
    sum = 0
    n_iters = length(num_arr)
    for i in 1:n_iters
        sum += num_arr[i]
    end
    return sum
end

f(num_arr)
```

`f(num_arr)` は配列(Array)以外にも適用可能だが `num_arr` という名前でいいのか？

```julia
f(1:5)
```

```julia
[typeof(1:5); [∘(fill(supertype, k)...)(typeof(1:5)) for k in 1:5]]
```

あと、`Base.sum` とローカル変数の名前が被っていても確かに問題ないが、それでもいいのか？

```julia
function f_rev1(num_abst_arr)
    total = 0
    n_iters = length(num_abst_arr)
    for i in 1:n_iters
        total += num_arr[i]
    end
    return total
end

f_rev1(num_arr)
```

`f_rev1([k for k in 1:5])` は計算できるが、`f_rev1(k for k in 1:5)` はエラーになる。

```julia
f_rev1([k for k in 1:5])
```

```julia
f_rev1(k for k in 1:5)
```

```julia
G = (k for k in 1:5)
```

```julia
function f_rev2(num_iter)
    total = 0
    for num in num_iter
        total += num
    end
    return total
end

f_rev2(num_arr)
```

```julia
f_rev2(k for k in 1:5)
```

引数 `num_iter` の要素が整数でない場合に型不安定！

```julia
@code_warntype f_rev2([1.0, 2.0, 3.0, 4.0, 5.0])
```

```julia
function f_rev3(num_iter)
    total = zero(first(num_iter))
    for num in num_iter
        total += num
    end
    return total
end

f_rev3(num_arr)
```

```julia
f_rev3(k for k in 1:5)
```

```julia
@code_warntype f_rev3([1.0, 2.0, 3.0, 4.0, 5.0])
```

`f_rev3(num_iter)` は引数 `num_iter` の要素が数でなくても使えるのに、`num_iter` という名前でいいのか？

```julia
f_rev3([[1, 2], [3, 4], [5, 6]])
```

```julia
function f_rev4(iter)
    total = zero(first(iter))
    for val in iter
        total += val
    end
    return total
end

f_rev4(num_arr)
```

```julia
f_rev4(k for k in 1:5)
```

```julia
f_rev4([[1, 2], [3, 4], [5, 6]])
```

Juliaでは全ての文が値を持つ式なので函数の終わりでの `return` は省略できる。

```julia
function f_rev5(iter)
    total = zero(first(iter))
    for val in iter
        total += val
    end
    total
end

f_rev5(num_arr)
```

```julia
f_rev5(k for k in 1:5)
```

```julia
f_rev5([[1, 2], [3, 4], [5, 6]])
```

`iter` が空の場合はエラーになる。

```julia
zeros(0)
```

```julia
f_rev5(zeros(0))
```

```julia
function f_rev6(iter, init = isempty(iter) ? zero(eltype(iter)) : zero(first(iter)))
    total = init
    for val in iter
        total += val
    end
    total
end

f_rev6(num_arr)
```

```julia
f_rev6(k for k in 1:5)
```

```julia
f_rev6([[1, 2], [3, 4], [5, 6]])
```

```julia
f_rev6(zeros())
```

```julia
f_rev6(Vector{Int64}[])
```

```julia
f_rev6(Vector{Int64}[], [0, 0])
```

函数の定義がこれだけ短いなら、「意味のある変数名」ではなく、1文字変数名でもいいんじゃね？

```julia
function f_rev7(A, o = isempty(A) ? zero(eltype(A)) : zero(first(A)))
    s = o
    for v in A
        s += v
    end
    s
end

f_rev7(num_arr)
```

```julia
f_rev7(k for k in 1:5)
```

```julia
f_rev7([[1, 2], [3, 4], [5, 6]])
```

```julia
f_rev7(zeros())
```

```julia
f_rev7(Vector{Int64}[], [0, 0])
```

```julia
A = [
    1 2
    3 4
]
f_rev7(A)
```

Juliaの配列のインデックスは任意始まりにできる。

```julia
using OffsetArrays
a = OffsetArray([k^3 for k in -3:3], -3:3)
```

```julia
@show collect(eachindex(a))

for i in eachindex(a)
    println("a[", i, "] = ", a[i])
end
```

```julia
@show pairs(a)
for (i, v) in pairs(a)
    println("a[", i, "] = ", v)
end
```

```julia
pairs(a)
```

```julia

```
