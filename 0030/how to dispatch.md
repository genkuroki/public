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
    display_name: Julia 1.7.2
    language: julia
    name: julia-1.7
---

```julia
module O

abstract type Mode end

struct Grid <: Mode end
read_grid(filename) = "read_grid: $filename"
readfile(filename, mode::Grid) = read_grid(filename)

DefaultMode() = Grid()
readfile(filename) = readfile(filename, DefaultMode())

struct PL3D <: Mode end
read_pl3d(filename) = "read_pl3d: $filename"
readfile(filename, mode::PL3D) = read_pl3d(filename)

struct Restart <: Mode end
read_restart(filename) = "read_restart: $filename"
readfile(filename, mode::Restart) = read_restart(filename)

struct Header <: Mode end
read_header(filename) = "read_header: $filename"
readfile(filename, mode::Header) = read_header(filename)

end
```

```julia
filename = "foo.ext"
@eval @show O.readfile($filename)
for mode in (O.Grid(), O.PL3D(), O.Restart(), O.Header())
    @eval @show O.readfile($filename, $mode)
end
```

上のようなスタイルでコードを書いておくと, モジュール O のコードを変更することなく, 他のモジュールで `O.readfile` メソッドを自由に拡張できるようになる.  このような拡張はJulia言語のエコシステムでは最も普通に行われていることである.  特に `Base.show` メソッドの独自拡張の頻度は極めて多い.

```julia
module P

using ..O

struct ModelP <: O.Mode end
read_modelp(filename) = "read_modelp: $filename"
O.readfile(filename, mode::ModelP) = read_modelp(filename)

end
```

```julia
filename = "foo.ext"
@eval @show O.readfile($filename)
for mode in (O.Grid(), O.PL3D(), O.Restart(), O.Header(), P.ModelP())
    @eval @show O.readfile($filename, $mode)
end
```

<!-- #region -->
## メモ1

```julia
readfile(filename::String, mode::Grid) = read_grid(filename)
```

と書くことはやめた方がよい. `filename::String` の部分に問題がある. このように書いてはいけない理由は文字列のように振る舞うが `String` 型でないオブジェクトが沢山あるからである. 例えば
<!-- #endregion -->

```julia
m = match(r"foo\d+.ext", "/bar/foo123.ext").match
```

```julia
my_readfile(filename::String, mode::O.Grid) = O.read_grid(filename)
my_readfile(m, O.Grid())
```

```julia
typeof(m)
```

```julia
m isa String
```

函数の引数の型を書くとすれば

```
readfile(filename::AbstractString, mode::Grid) = read_grid(filename)
```

と `String` ではなく `AbstractString` を使うべき.

```julia
m isa AbstractString
```

```julia
my_readfile(filename::AbstractString, mode::O.Grid) = O.read_grid(filename)
my_readfile(m, O.Grid())
```

<!-- #region -->
## メモ2

さらに module O の中で

```julia
readfile(filename, mode::Mode) = println("Not implemented for ", typeof(mode))
```

のように定義すると, デバッグが困難になる原因になるので注意が必要. 詳しくは

* https://www.oxinabox.net/2020/04/19/Julia-Antipatterns.html#notimplemented-exceptions

を参照.

たとえば, 函数の引数の型注釈を

```
readfile(filename::String, mode::Grid) = read_grid(filename)
```

と書いているときに(これはメモ1で述べたようにやめた方がよい), 上のように not implemented エラーを出すようにしてあると, 以下のようなことが起こって困惑してしまうことになる.
<!-- #endregion -->

```julia
module O1

abstract type Mode end
readfile(filename, mode::Mode) = println("Not implemented for ", typeof(mode))

struct Grid <: Mode end
read_grid(filename) = "read_grid: $filename"
readfile(filename::String, mode::Grid) = read_grid(filename)

end
```

```julia
O1.readfile(m, O1.Grid())
```

```julia
m
```

```julia
typeof(m)
```

以下のように Not implemented エラーのために函数を定義していない方がエラーメッセージが分かり易い.

```julia
module O2

abstract type Mode end

struct Grid <: Mode end
read_grid(filename) = "read_grid: $filename"
readfile(filename::String, mode::Grid) = read_grid(filename)

end
```

```julia
@show m
O2.readfile(m, O2.Grid())
```

```julia

```
