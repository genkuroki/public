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
    display_name: Julia current stable release
    language: julia
    name: julia
---

# ColabでJulia言語を使った統計学の勉強の仕方

* 黒木玄
* 2025-05-13, 2025-05-21
$
\newcommand\op{\operatorname}
\newcommand\R{{\mathbb R}}
\newcommand\Z{{\mathbb Z}}
\newcommand\var{\op{var}}
\newcommand\std{\op{std}}
\newcommand\eps{\varepsilon}
\newcommand\T[1]{T_{(#1)}}
\newcommand\bk{\bar\kappa}
\newcommand\X{{\mathscr X}}
$

このノートブックは[Google Colabで実行できる](https://colab.research.google.com/github/genkuroki/public/blob/main/0054/How%20to%20use%20Julia%20language%20in%20Google%20Colab%20for%20learning%20statistics.ipynb).

<!-- #region toc=true -->
<h1>目次<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#Google-ColabでのJulia言語の使い方" data-toc-modified-id="Google-ColabでのJulia言語の使い方-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Google ColabでのJulia言語の使い方</a></span><ul class="toc-item"><li><span><a href="#ColabでのJuliaの実行" data-toc-modified-id="ColabでのJuliaの実行-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>ColabでのJuliaの実行</a></span></li><li><span><a href="#グラフの描き方" data-toc-modified-id="グラフの描き方-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>グラフの描き方</a></span></li><li><span><a href="#標準正規分布乱数のプロット" data-toc-modified-id="標準正規分布乱数のプロット-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>標準正規分布乱数のプロット</a></span></li><li><span><a href="#確率分布の扱い方" data-toc-modified-id="確率分布の扱い方-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>確率分布の扱い方</a></span></li><li><span><a href="#正規分布の確率密度関数のプロット" data-toc-modified-id="正規分布の確率密度関数のプロット-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>正規分布の確率密度関数のプロット</a></span></li></ul></li><li><span><a href="#Anscombeの例のプロット" data-toc-modified-id="Anscombeの例のプロット-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Anscombeの例のプロット</a></span><ul class="toc-item"><li><span><a href="#RDatasets.jlパッケージのインストール" data-toc-modified-id="RDatasets.jlパッケージのインストール-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>RDatasets.jlパッケージのインストール</a></span></li><li><span><a href="#データのプロットの仕方" data-toc-modified-id="データのプロットの仕方-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>データのプロットの仕方</a></span></li></ul></li><li><span><a href="#Datasaurusの散布図のプロット" data-toc-modified-id="Datasaurusの散布図のプロット-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Datasaurusの散布図のプロット</a></span><ul class="toc-item"><li><span><a href="#データの取得" data-toc-modified-id="データの取得-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>データの取得</a></span></li><li><span><a href="#散布図の作成" data-toc-modified-id="散布図の作成-3.2"><span class="toc-item-num">3.2&nbsp;&nbsp;</span>散布図の作成</a></span></li></ul></li><li><span><a href="#中心極限定理のプロット" data-toc-modified-id="中心極限定理のプロット-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>中心極限定理のプロット</a></span><ul class="toc-item"><li><span><a href="#素朴なワークフロー" data-toc-modified-id="素朴なワークフロー-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>素朴なワークフロー</a></span></li><li><span><a href="#関数を作ろう" data-toc-modified-id="関数を作ろう-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>関数を作ろう</a></span></li><li><span><a href="#問題:-自分で関数を定義して実行してみよ." data-toc-modified-id="問題:-自分で関数を定義して実行してみよ.-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>問題: 自分で関数を定義して実行してみよ.</a></span></li></ul></li><li><span><a href="#統計モデルの例" data-toc-modified-id="統計モデルの例-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>統計モデルの例</a></span><ul class="toc-item"><li><span><a href="#t分布が得られる統計モデルの例" data-toc-modified-id="t分布が得られる統計モデルの例-5.1"><span class="toc-item-num">5.1&nbsp;&nbsp;</span>t分布が得られる統計モデルの例</a></span></li><li><span><a href="#t分布が得られる統計モデルの例の別バージョン" data-toc-modified-id="t分布が得られる統計モデルの例の別バージョン-5.2"><span class="toc-item-num">5.2&nbsp;&nbsp;</span>t分布が得られる統計モデルの例の別バージョン</a></span></li><li><span><a href="#負の二項分布が得られる統計モデルの例" data-toc-modified-id="負の二項分布が得られる統計モデルの例-5.3"><span class="toc-item-num">5.3&nbsp;&nbsp;</span>負の二項分布が得られる統計モデルの例</a></span></li></ul></li></ul></div>
<!-- #endregion -->

```julia
# Google Colabと自分のパソコンの両方で使えるようにするための工夫

import Pkg

"""すでにPkg.add済みのパッケージのリスト (高速化のために用意)"""
_packages_added = [info.name for (uuid, info) in Pkg.dependencies() if info.is_direct_dep]

"""_packages_added内にないパッケージをPkg.addする"""
add_pkg_if_not_added_yet(pkg) = if !(pkg in _packages_added)
    println(stderr, "# $(pkg).jl is not added yet, so let's add it.")
    Pkg.add(pkg)
end

"""expr::Exprからusing内の`.`を含まないモジュール名を抽出"""
function find_using_pkgs(expr::Expr)
    pkgs = String[]
    function traverse(expr::Expr)
        if expr.head == :using
            for arg in expr.args
                if arg.head == :. && length(arg.args) == 1
                    push!(pkgs, string(arg.args[1]))
                elseif arg.head == :(:) && length(arg.args[1].args) == 1
                    push!(pkgs, string(arg.args[1].args[1]))
                end
            end
        else
            for arg in expr.args arg isa Expr && traverse(arg) end
        end
    end
    traverse(expr)
    pkgs
end

"""必要そうなPkg.addを追加するマクロ"""
macro autoadd(expr)
    pkgs = find_using_pkgs(expr)
    :(add_pkg_if_not_added_yet.($(pkgs)); $expr)
end

using Random
Random.seed!(4649373)

@autoadd begin
using Distributions
using RDatasets
using StatsPlots
default(fmt=:png)
end
```

<!-- #region -->
上のセルでは Distributions.jl, RDataset.jl, StatsPlots.jl の3つのパッケージをインストール(`Pkg.add`)した後に, それらの `using` を実行している. 必要な `Pkg.add` のコードは `@autoadd` マクロが追加してくれている.

存在するパッケージ A.jl について `using A` を実行したとき,

```
ArgumentError: Package A not found in current path.
- Run `import Pkg; Pkg.add("A")` to install the A package.
```

と表示されたならば

```julia
import Pkg; Pkg.add("A")
using A
```

または, `@autoadd` マクロを使って

```julia
@autoadd using A
```

を実行すればよい.

<a href="https://docs.julialang.org/en/v1/">Julia言語</a>のBaseやStandard Libraryに含まれるパッケージは`Pkg.add`をしなくても使用できる. たとえば上のセルのように `using Random` には `@autoadd` を適用する必要はない. 他にも `using LinearAlgebra` や `using Printf` の類も `@autoadd` は必要ない. 
<!-- #endregion -->

## Google ColabでのJulia言語の使い方

<!-- #region -->
### ColabでのJuliaの実行

(1) ブラウザでGoogleアカウントのどれかにログインしておきます.

(2) [Google Colab](https://colab.research.google.com/)にアクセスする.

(3) 「ノートブックを開く」の「GitHub」を選択する.

(4) GitHubにおいてある `ipynb` ファイルのURLを入力してEnterキーを押す.  例えば

 * `https://github.com/genkuroki/public/blob/main/0054/How%20to%20use%20Julia%20language%20in%20Google%20Colab%20for%20learning%20statistics.ipynb`

というURLを入力する.

(5) 実際にその例のURLを入力してEnterキーを押すと, このファイルがGoogle Colabで開かれる.

(6) そのノートブックの全体をColabで実行し直したければ, 「ランタイム」→「すべてのセルを実行」を選択する.

(7) 適当にGoogle Colabの使い方を検索して調べればより詳しい使い方が分かる.


* 各セルの先頭に `?` と入力した後に関数名などを入れるとヘルプを読むことができる.
* 各セルの先頭に `]` と入力した後にパッケージ管理モードのコマンドを入力して実行できる.
* タブキーによる補完を使える.
* 各セルの最後に `;` を付けて実行すると計算結果が表示されない.

__問題:__ 以上を実際に行ってみよ.
<!-- #endregion -->

```julia
1 + 1
```

```julia
sin(pi/6)
```

```julia
sinpi(1/6)
```

<!-- #region -->
### グラフの描き方

(8) Colabで統計学対応のグラフ作画パッケージを使うためには次を実行する:

```julia
import Pkg
Pkg.add("StatsPlots")
using StatsPlots
```

このノートブックでは最初のセルでこれと同等のことを実行できるようにしてあるので, 最初のセルを実行しておけばよい.
<!-- #endregion -->

```julia
plot(sin)
```

```julia
plot(Normal(3, 10))
```

### 標準正規分布乱数のプロット

```julia
# 標準正規分布の乱数を10^4個生成
Z = randn(10^4);
```

```julia
histogram(Z; norm=true, alpha=0.5, label="")
```

```julia
plot!(x -> exp(-x^2/2)/sqrt(2pi), -4, 4; label="", lw=3)
```

<!-- #region -->
### 確率分布の扱い方

(9) 確率分布を扱うためのパッケージを使うためには次を実行する:

```julia
import Pkg
Pkg.add("StatsPlots")
using StatsPlots
```

このノートブックでは最初のセルでこれと同等のことを実行できるようにしてあるので, 最初のセルを実行しておけばよい.
<!-- #endregion -->

```julia
dist = Binomial(20, 0.3)
bar(dist; alpha=0.5, label="Binomial(20, 0.3)")
```

### 正規分布の確率密度関数のプロット

```julia
X = rand(Normal(2, 3), 10^4);
```

```julia
histogram(X; norm=true, alpha=0.5, label="sample of Normal(2, 3)")
```

```julia
plot!(Normal(2, 3); label="Normal(2, 3)", lw=2)
```

## Anscombeの例のプロット

<!-- #region -->
### RDatasets.jlパッケージのインストール

Rのデータセット達を扱うためのパッケージを入れるためには次を実行する:

```julia
import Pkg
Pkg.add("RDatasets")
using RDatasets
```

このノートブックでは最初のセルでこれと同等のことを実行できるようにしてあるので, 最初のセルを実行しておけばよい.
<!-- #endregion -->

```julia
anscombe = dataset("datasets", "anscombe")
```

### データのプロットの仕方

以下ではデータ1の場合のプロットの仕方を説明しよう.

```julia
# x, y にデータを入れる
x, y = anscombe.X1, anscombe.Y1
```

```julia
# 散布図を描いてみる
using StatsPlots
scatter(x, y)
```

```julia
# xlim, ylimなどを追加
scatter(x, y; label="", xlim=(3, 20), ylim=(2, 14))
```

```julia
# データの標本平均や不偏分散・不偏共分散を計算
xbar = mean(x)
```

```julia
ybar = mean(y)
```

```julia
sx2 = var(x)
```

```julia
sy2 = var(y)
```

```julia
sxy = cov(x, y)
```

```julia
betahat = sxy/sx2
```

```julia
alphahat = ybar - betahat*xbar
```

```julia
scatter(x, y; label="", xlim=(3, 20), ylim=(2, 14))
plot!(x -> alphahat + betahat*x, 3, 20; label="", lw=2)
```

```julia
scatter(x, y; label="", xlim=(3, 20), ylim=(2, 14), title="Anscombe 1")
plot!(x -> alphahat + betahat*x, 3, 20; label="", lw=2, ls=:dash)
```

```julia
# design matrix
X = x .^ (0:1)'
```

```julia
# 最小二乗法を一発実現 (計画行列の一般逆行列をyにかける)
alphahat2, betahat2 = X \ y
```

```julia
# 2つの直線はぴったり重なり合う.
scatter(x, y; label="", xlim=(3, 20), ylim=(2, 14), title="Anscombe 1")
plot!(x -> alphahat + betahat*x, 3, 20; label="", lw=2)
plot!(x -> alphahat2 + betahat2*x, 3, 20; label="", lw=2, ls=:dash)
```

__問題:__ 他のアンスコムのデータについて同様のグラフを作成せよ.


## Datasaurusの散布図のプロット

以下のデータは「[条件付き確率分布, 尤度, 推定, 記述統計](https://nbviewer.org/github/genkuroki/Statistics/blob/master/2022/06%20Conditional%20distribution%2C%20likelihood%2C%20estimation%2C%20and%20summary.ipynb)」からのコピー＆ペースト.


### データの取得

<!--
* http://www.thefunctionalart.com/2016/08/download-datasaurus-never-trust-summary.html
  * https://www.dropbox.com/sh/xaxpz3pm5r5awes/AADUbGVagF9i4RmM9JkPtviEa?dl=0
-->
* https://www.dropbox.com/sh/xaxpz3pm5r5awes/AADUbGVagF9i4RmM9JkPtviEa?dl=0
* https://visualizing.jp/the-datasaurus-dozen/
* https://www.openintro.org/data/index.php?data=datasaurus

```julia
datasaurus = [
    55.3846 97.1795
    51.5385 96.0256
    46.1538 94.4872
    42.8205 91.4103
    40.7692 88.3333
    38.7179 84.8718
    35.6410 79.8718
    33.0769 77.5641
    28.9744 74.4872
    26.1538 71.4103
    23.0769 66.4103
    22.3077 61.7949
    22.3077 57.1795
    23.3333 52.9487
    25.8974 51.0256
    29.4872 51.0256
    32.8205 51.0256
    35.3846 51.4103
    40.2564 51.4103
    44.1026 52.9487
    46.6667 54.1026
    50.0000 55.2564
    53.0769 55.6410
    56.6667 56.0256
    59.2308 57.9487
    61.2821 62.1795
    61.5385 66.4103
    61.7949 69.1026
    57.4359 55.2564
    54.8718 49.8718
    52.5641 46.0256
    48.2051 38.3333
    49.4872 42.1795
    51.0256 44.1026
    45.3846 36.4103
    42.8205 32.5641
    38.7179 31.4103
    35.1282 30.2564
    32.5641 32.1795
    30.0000 36.7949
    33.5897 41.4103
    36.6667 45.6410
    38.2051 49.1026
    29.7436 36.0256
    29.7436 32.1795
    30.0000 29.1026
    32.0513 26.7949
    35.8974 25.2564
    41.0256 25.2564
    44.1026 25.6410
    47.1795 28.7180
    49.4872 31.4103
    51.5385 34.8718
    53.5897 37.5641
    55.1282 40.6410
    56.6667 42.1795
    59.2308 44.4872
    62.3077 46.0256
    64.8718 46.7949
    67.9487 47.9487
    70.5128 53.7180
    71.5385 60.6410
    71.5385 64.4872
    69.4872 69.4872
    46.9231 79.8718
    48.2051 84.1026
    50.0000 85.2564
    53.0769 85.2564
    55.3846 86.0256
    56.6667 86.0256
    56.1538 82.9487
    53.8462 80.6410
    51.2821 78.7180
    50.0000 78.7180
    47.9487 77.5641
    29.7436 59.8718
    29.7436 62.1795
    31.2821 62.5641
    57.9487 99.4872
    61.7949 99.1026
    64.8718 97.5641
    68.4615 94.1026
    70.7692 91.0256
    72.0513 86.4103
    73.8462 83.3333
    75.1282 79.1026
    76.6667 75.2564
    77.6923 71.4103
    79.7436 66.7949
    81.7949 60.2564
    83.3333 55.2564
    85.1282 51.4103
    86.4103 47.5641
    87.9487 46.0256
    89.4872 42.5641
    93.3333 39.8718
    95.3846 36.7949
    98.2051 33.7180
    56.6667 40.6410
    59.2308 38.3333
    60.7692 33.7180
    63.0769 29.1026
    64.1026 25.2564
    64.3590 24.1026
    74.3590 22.9487
    71.2821 22.9487
    67.9487 22.1795
    65.8974 20.2564
    63.0769 19.1026
    61.2821 19.1026
    58.7179 18.3333
    55.1282 18.3333
    52.3077 18.3333
    49.7436 17.5641
    47.4359 16.0256
    44.8718 13.7180
    48.7179 14.8718
    51.2821 14.8718
    54.1026 14.8718
    56.1538 14.1026
    52.0513 12.5641
    48.7179 11.0256
    47.1795  9.8718
    46.1538  6.0256
    50.5128  9.4872
    53.8462 10.2564
    57.4359 10.2564
    60.0000 10.6410
    64.1026 10.6410
    66.9231 10.6410
    71.2821 10.6410
    74.3590 10.6410
    78.2051 10.6410
    67.9487  8.7180
    68.4615  5.2564
    68.2051  2.9487
    37.6923 25.7692
    39.4872 25.3846
    91.2821 41.5385
    50.0000 95.7692
    47.9487 95.0000
    44.1026 92.6923
];
```

### 散布図の作成

```julia
# 行列Aの第j列はA[:,j]
@show datasaurus[:,1];
```

```julia
scatter(datasaurus[:,1], datasaurus[:,2]; label="", title="Datasaurus")
```

__問題:__ [Datasaurusについて検索](https://www.google.com/search?q=Datasaurus)して見つけた解説を読め.


## 中心極限定理のプロット


### 素朴なワークフロー

以下のセルの内容を julia の `julia> ` プロンプトに順番に入力すれば(コピー＆ペーストすれば)同じ結果が得られる.  各行の最後にセミコロン `;` を追加すれば計算結果の出力を抑制できる.

```julia
# 確率分布を dist と書く.
dist = Gamma(2, 3)
```

```julia
# 確率分布 dist のサイズ n のサンプルを L 個生成
n = 10
L = 10^4
Xs = [rand(dist, n) for _ in 1:L]
```

```julia
# L個のサイズnのサンプルの各々の標本平均を計算
Xbars = mean.(Xs)
```

```julia
# Xbarのヒストグラムを表示
histogram(Xbars; norm=true, alpha=0.5, label="", title="$dist, n=$n")
```

```julia
# 中心極限定理による正規分布近似を設定
mu = mean(dist)
sigma = std(dist)
normal_approx = Normal(mu, sigma/sqrt(n))
```

```julia
# 上のグラフに重ねて正規分布をプロット
plot!(normal_approx; label="normal approx", lw=2)
```

$n=10$ が小さすぎてずれが大きい.

```julia
# nを大きくしてやり直してみる.
n = 100
L = 10^4
Xs = [rand(dist, n) for _ in 1:L]
Xbars = mean.(Xs)
histogram(Xbars; norm=true, alpha=0.5, label="", title="$dist, n=$n")
mu = mean(dist)
sigma = std(dist)
normal_approx = Normal(mu, sigma/sqrt(n))
plot!(normal_approx; label="normal approx", lw=2)
```

$n=100$ にしたら, 正規分布とよく一致するようになった.

```julia
# Lも大きくしてやり直してみる.
n = 100
L = 10^5
Xs = [rand(dist, n) for _ in 1:L]
Xbars = mean.(Xs)
histogram(Xbars; norm=true, alpha=0.5, label="", title="$dist, n=$n")
mu = mean(dist)
sigma = std(dist)
normal_approx = Normal(mu, sigma/sqrt(n))
plot!(normal_approx; label="normal approx", lw=2)
```

```julia
# Lも大きくしてやり直してみる.
n = 1000
L = 10^5
Xs = [rand(dist, n) for _ in 1:L]
Xbars = mean.(Xs)
histogram(Xbars; norm=true, alpha=0.5, label="", title="$dist, n=$n")
mu = mean(dist)
sigma = std(dist)
normal_approx = Normal(mu, sigma/sqrt(n))
plot!(normal_approx; label="normal approx", lw=2)
```

### 関数を作ろう

上のように素朴に毎回コードを入力することは非常に面倒である.

同じような仕事を繰り返したい場合には, その仕事を関数化して1行の入力の繰り返しで実行できるようにしておく方がよい.


### 問題: 自分で関数を定義して実行してみよ.

以下のセルのように関数を定義しておくと, 同じような仕事を何度も楽に実行できるようになる.

```julia
using Distributions
using StatsPlots
default(fmt=:png)

function hello_sine()
    println("Hello, Sine!")
    plot(sin; label="y=sin(x)")
end

mypdf(dist, x) = pdf(dist, x)
mypdf(dist::ContinuousUnivariateDistribution, x) =
    minimum(dist) < x < maximum(dist) ? pdf(dist, x) : zero(eltype(dist))
mypdf(dist::DiscreteUnivariateDistribution, x) = pdf(dist, round(Int, x))
mydistname(dist) = replace(string(dist), r"{[^}]*}"=>"")
mydistname(dist::InverseGamma) = "InverseGamma(α=$(shape(dist)), θ=$(scale(dist)))"

function plot_dist(dist; xlim0=nothing)
    distname = mydistname(dist)
    if isnothing(xlim0)
        mu = mean(dist)
        sigma = std(dist)
        a = max(minimum(dist), mu - 4.5sigma)
        b = min(maximum(dist), mu + 4.5sigma)
        if dist isa DiscreteUnivariateDistribution
            a, b = round(a)-1, round(b)+1
        else
            a, b = a-0.025(b-a), b+0.025(b-a)
        end
        xlim0 = (a, b)
    end
    xs = range(xlim0..., 1000)
    plot(xs, x -> mypdf(dist, x); label="", title="$distname")
    plot!(; size=(400, 250))
end

function plot_central_limit_theorem(dist, n; L=10^5, bin=:auto, xlim1=nothing)
    distname = mydistname(dist)
    mu = mean(dist)
    sigma = std(dist)
    Xs = [rand(dist, n) for _ in 1:L]
    Xbars = mean.(Xs)
    normal_approx = Normal(mu, sigma/sqrt(n))
    
    if dist isa DiscreteUnivariateDistribution
        mu = mean(dist)
        sigma = std(dist)
        a = round(n*mu - 4.5sqrt(n)*sigma)
        b = round(n*mu + 4.5sqrt(n)*sigma)
        ran = a-0.5:b+0.5
        bin = ran / n
    end
    
    #histogram(Xbars; bin, norm=true, alpha=0.5, label="Xbars")
    stephist(Xbars; bin, norm=true, label="Xbars")
    plot!(normal_approx; lw=2, label="normal approx")
    isnothing(xlim1) || plot!(; xlim=xlim1)
    title!("$distname, n=$n")
    plot!(; size=(400, 250))
end

function plot_dist_clt(dist, n; L=10^5, xlim0=nothing, xlim1=nothing, bin=100)
    P0 = plot_dist(dist; xlim0)
    P1 = plot_central_limit_theorem(dist, n; L, bin, xlim1)
    plot(P0, P1; size=(800, 250), layout=(1, 2))
    plot!(; titlefontsize=10)
end
```

```julia
hello_sine()
```

以下で登場する確率分布(`Uniform()`, etc.)の解説が以下の場所にある:

* https://juliastats.org/Distributions.jl/stable/univariate/#Continuous-Distributions
* https://juliastats.org/Distributions.jl/stable/univariate/#Discrete-Distributions

```julia
plot_dist_clt(Uniform(), 10)
```

```julia
plot_dist_clt(Laplace(), 10)
```

```julia
plot_dist_clt(Gamma(2, 3), 10)
```

```julia
plot_dist_clt(Gamma(2, 3), 100)
```

```julia
plot_dist_clt(Gamma(2, 3), 1000)
```

```julia
plot_dist_clt(Exponential(), 10)
```

```julia
plot_dist_clt(Exponential(), 100)
```

```julia
plot_dist_clt(LogNormal(), 10)
```

```julia
plot_dist_clt(LogNormal(), 100)
```

```julia
plot_dist_clt(InverseGamma(3, 3), 10; xlim0=(-0.2, 8.2), xlim1=(-0.5, 3.5))
```

```julia
plot_dist_clt(InverseGamma(3, 3), 100; xlim0=(-0.2, 6.2), xlim1=(0.9, 2.1))
```

```julia
plot_dist_clt(InverseGamma(3, 3), 1000; xlim0=(-0.2, 6.2), xlim1=(1.3, 1.7))
```

```julia
plot_dist_clt(Bernoulli(), 10)
```

```julia
plot_dist_clt(Bernoulli(), 100)
```

```julia
plot_dist_clt(Bernoulli(0.2), 10)
```

```julia
plot_dist_clt(Bernoulli(0.2), 100)
```

```julia
plot_dist_clt(Poisson(), 10)
```

```julia
plot_dist_clt(Poisson(), 100)
```

## 統計モデルの例

階層的な統計モデルは「乱数生成の連鎖」というイメージで理解しておくと理解し易い.


### t分布が得られる統計モデルの例

たとえば, 

$$
\begin{aligned}
& Z \sim \op{Normal}(0, 1) \\
& W \sim \op{Chisq}(\nu) \\
& T = \frac{Z}{\sqrt{W/\nu}}
\end{aligned}
$$

は以下のような統計モデルを表す:

1. 標準正規分布に従う乱数 $Z$ を生成する.
2. 自由度 $\nu$ のχ²分布に従う乱数 $W$ を生成する.
3. $T = Z\big/\sqrt{W/\nu}$ と定める.

このとき, $T$ は自由度 $\nu$ のt分布に従う乱数になる.

これを確認しよう.

```julia
using Distributions
using StatsPlots
default(fmt=:png)

function plot_t(T, nu; Tlabel="T")
    stephist(T; norm=true, label=Tlabel)
    plot!(TDist(nu); label="TDist(nu)", ls=:dash)
    plot!(xlim=(-6, 6))
    title!("nu = $nu")
end

function plot_model_t1(; nu=3, L=10^6)
    Z = rand(Normal(0, 1), L)
    W = rand(Chisq(nu), L)
    T = @. Z / sqrt(W / nu)
    plot_t(T, nu; Tlabel="T = Z/√(W/ν)")
    plot!(; size=(500, 300), titlefontsize=14, legendfontsize=10)
end
```

```julia
plot_model_t1(; nu=3)
```

```julia
plot_model_t1(; nu=10)
```

### t分布が得られる統計モデルの例の別バージョン

$$
\begin{aligned}
& Z \sim \op{Normal}(0, 1) \\
& W \sim \op{Chisq}(\nu) \\
& T = \frac{Z}{\sqrt{W/\nu}}
\end{aligned}
$$

というモデルにおいて, $T$ は期待値 $0$ と分散 $1/(W/\nu)$ を持つ正規分布の乱数であるとも解釈できる. 
この解釈によって次の統計モデルを考える.

$$
\begin{aligned}
& W \sim \op{Chisq}(\nu) \\
& T \sim \op{Normal}\left(0, 1\big/\sqrt{W/\nu}\right)
\end{aligned}
$$

この $T$ も自由度 $\nu$ のt分布に従う乱数になる.

これを確認しよう.

```julia
function plot_model_t2(; nu=3, L=10^6)
    W = rand(Chisq(nu), L)
    T = @. rand(Normal(0, 1/sqrt(W/nu)))
    plot_t(T, nu; Tlabel="T~N(0, 1/√(W/ν))")
end
```

```julia
plot_model_t2(; nu=3)
```

```julia
plot_model_t2(; nu=10)
```

### 負の二項分布が得られる統計モデルの例

$$
\begin{aligned}
& \Lambda \sim \op{Gamma}(\alpha, \theta) \\
& K \sim \op{Poisson}(\Lambda)
\end{aligned}
$$

このとき $K \sim \op{NegativeBinomial}(\alpha, 1/(1+\theta))$ となることを確認しよう.

```julia
function plot_gampoi(α, θ; L = 10^6)
    p = 1/(1 + θ)
    Λ = rand(Gamma(α, θ), L) # ガンマ分布で乱数達を大量生成
    M_gampoi = @. rand(Poisson(Λ)) # その各々を期待値とするPoisson分布の乱数を生成
    M_negbin = rand(NegativeBinomial(α, p), L) # 直接的に負の二項分布の乱数を大量生成

    # 比較のための同時プロット
    binmin, binmax = round.(quantile.(Ref(M_negbin), (0.001, 0.999)))
    stephist(M_gampoi; norm=true, bin=binmin-0.5:binmax+0.5, label="Gamma-Poisson")
    stephist!(M_negbin; norm=true, bin=binmin-0.5:binmax+0.5, ls=:dash, label="NegativeBinomial")
    title!("α=$α, θ=$θ")
end
```

```julia
plot_gampoi(2, 1)
```

```julia
plot_gampoi(2, 3)
```

```julia
plot_gampoi(2, 10)
```

```julia

```
