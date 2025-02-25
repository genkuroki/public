# Julia言語のPlots.jl pgfplotsx()で日本語を使う方法

* 黒木玄
* 2025-02-12

## ファイルの説明

`pgfplotsx-ja-example.ipynb` はグラフが含まれる `pgfplotsx-ja-bxjsarticle.pdf`, `pgfplotsx-ja-bxjsarticle.tex`, `pgfplotsx-ja-standalone.tex` を作成するためのJuliaカーネルのJupyter notebookである.

`pgfplotsx-ja-bxjsarticle.tex` を LuaLaTeX でコンパイルすると, `pgfplotsx-ja-bxjsarticle.pdf`, `pgfplotsx-ja-bxjsarticle.tex` が読み込まれ, 2つの図を含む `pgfplotsx-ja-bxjsarticle.pdf` が出力される.

`pgfplotsx-ja-standalone.tex` を LuaLaTeX でコンパイルすると, 図のPDFファイル `pgfplotsx-ja-standalone.pdf` が出力される.

## Plots.jl pgfplotsx()で日本語を使うためのポイント

```julia
using Plots
pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE=[raw"\usepackage{luatexja}"]
@eval Plots pgfx_sanitize_string(s::AbstractString) = s
```

のように設定すればよい. 

日本語を使うために `PGFPlotsX.CUSTOM_PREAMBLE=[raw"\usepackage{luatexja}"]` が必要なことはすぐにわかる.

非自明なのは, `Plots.pgfx_sanitize_string(s::AbstractString)` を無効にすることである. 

`Plots.pgfx_sanitize_string(s::AbstractString)`を無効にしておかないと, その副作用が大き過ぎて title, legend, guide などでうまいこと日本語を使えなくなってしまう.

`Plots.pgfx_sanitize_string(s::AbstractString)`を無効にすることによって, title, legend, guide などでLaTeXのコードを使えるようになる.

実際に動くコードについては `pgfplotsx-ja-example.ipynb` を参照せよ.

## 非破壊的な方法

```julia
using LaTeXStrings
using Plots
pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE=[raw"\usepackage{luatexja}"]
```

として, 文字列を`LaTeXString`でラップするようにすれば, `pgfx_sanitize_string(s::AbstractString)`を無効にする破壊的な変更なしで, pgfplotsx()で日本語を使える.

この方法も `pgfplotsx-ja-example.ipynb` で紹介してある.
