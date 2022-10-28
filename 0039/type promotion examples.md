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
    display_name: Julia 1.8.2
    language: julia
    name: julia-1.8
---

https://twitter.com/taketo1024_2/status/1585942690495336449

```julia
using LinearAlgebra: LinearAlgebra
using SparseArrays
using Polynomials
```

```julia
Base.promote_op(LinearAlgebra.matprod, Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h})
```

これ↑が以下のようになる理由.

```julia
p = Polynomial(Rational[0, 1], :h)
```

```julia
eltype(p)
```

```julia
typeof(p)
```

`Rational` は抽象型なので `Polynomial{Rational, :h}` 型のオブジェクトを使用すると型不安定になって計算速度が劣化する危険性があるので要注意.  Juliaでのプログラミングでは具象型の情報が正しく伝搬するように書くことが最も重要である.  抽象型と具象型の区別は非常に重要. 

抽象型を避けるべき場合があることについては公式ドキュメントの

* https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-abstract-container

を読めばわかる.

Juliaの大きな特徴の1つは, 多重ディスパッチを使ったコードによって, 数値の型が混合している場合の演算を定義していることである.  この [mixed type arithmetic 実装](https://www.google.com/search?q=mixed-type+site%3Adocs.julialang.org)は多重ディスパッチを採用していない場合にはかなりやっかいな問題になる.

そのための型のプロモーションのルールはJuliaコードできちんと記述されているので, Juliaの型の取り扱いは実際には全然「ゆるふわ」ではない.  むしろ「ゆるふわ」でないせいでJulia初心者は失敗することが多いように思われる.

```julia
a = sparse([p;;])
```

```julia
-a
```

`-a` の `eltype` が `Polynomial{Rational, :h}` ではなく具象型の `Polynomial{Rational{Int64}, :h}` になっている.

```julia
(-a)*a
```

`(-a)*a` の `eltype` が `Polynomial{Any, :h}` になってしまった!


こうなってしまった理由はすでに説明したように次のように型がプロモートするようになっているからである.

```julia
Base.promote_op(LinearAlgebra.matprod, Polynomial{Rational{Int}, :h}, Polynomial{Rational, :h})
```

ただし, `Any` が出て来るより前に `Rational` という抽象型を `eltype` とする多項式を使ってしまっているので, 計算速度的にはすでに損失が生じていたと考えられる.


こうならないようにするためには, `p = Polynomial(Rational[0, 1], :h)` を次に置き換えればよい(`Rational` → `Rational{Int}`).

```julia
q = Polynomial(Rational{Int}[0, 1], :h)
```

```julia
typeof(q)
```

```julia
eltype(q)
```

```julia
isconcretetype(Rational) # Rationalは具象型ではない
```

```julia
isconcretetype(Rational{Int}) # Rational{Int}は具象型である
```

```julia
b = sparse([q;;])
```

```julia
-b
```

```julia
(-b)*b
```

```julia
Base.promote_op(LinearAlgebra.matprod, Polynomial{Rational{Int}, :h}, Polynomial{Rational{Int}, :h})
```

他の例も示しておこう.

```julia
r = Polynomial(Rational{BigInt}[0, 1], :h)
```

```julia
typeof(r)
```

```julia
c = sparse([r;;])
```

```julia
(-b)*c
```

```julia
Base.promote_op(LinearAlgebra.matprod, Polynomial{Rational{Int}, :h}, Polynomial{Rational{BigInt}, :h})
```

このように `Int` と `BigInt` の演算では `BigInt` の結果が得られるようになっている.

```julia

```
