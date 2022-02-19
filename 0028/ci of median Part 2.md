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

# 中央値の信頼区間2

Part 1: [中央値の信頼区間](https://github.com/genkuroki/public/blob/main/0028/confidence%20interval%20of%20median.ipynb)
$
\newcommand\R{\Bbb{R}}
\newcommand\Beta{\operatorname{Beta}}
\newcommand\Binomial{\operatorname{Binomial}}
\newcommand\dist{\operatorname{dist}}
\newcommand\empirical{\operatorname{empirical}}
\newcommand\pval{\operatorname{pval}}
\newcommand\median{\operatorname{median}}
\newcommand\ci{\operatorname{ci}}
\newcommand\pdf{\operatorname{pdf}}
\newcommand\cdf{\operatorname{cdf}}
\newcommand\ecdf{\operatorname{ecdf}}
\newcommand\quantile{\operatorname{quantile}}
\newcommand\on{\operatorname}
$

```julia
ENV["COLUMNS"] = 200

using Distributions
using StatsPlots
default(titlefontsize=10, fmt=:png)
using Random
using StatsBase
using QuadGK
using StaticArrays
using DataFrames

name(dist::UnivariateDistribution) = replace(string(dist), r"{[^{.]*}"=>"")
```

```julia
function pdf_median_true(n, z)
    0 < z < 1 || return 0.0
    m = n / 2
    p(x, y) = pdf(Dirichlet(SVector(m, 1, m)), SVector(x, y-x, 1-y))
    2quadgk(x -> p(x, 2z - x), 0, z)[1]
end

function plot_mediandists(n; kwargs...)
    plot(x -> pdf_median_true(n, x), 0, 1; label="true dist")
    plot!(x -> pdf(Beta((n+1)/2, (n+1)/2), x), 0, 1; ls=:dash, label="n")
    plot!(x -> pdf(Beta((n+2)/2, (n+2)/2), x), 0, 1; ls=:dashdot, label="n'=n+1")
    title!("dist of median and approx: n = $n")
    plot!(; kwargs...)
end

function plot_mediandist_approx()
    PP = []
    P = plot_mediandists(2; legend=:bottom)
    push!(PP, P)
    for n in 4:2:18
    P = plot_mediandists(n; legend=false)
        push!(PP, P)
    end
    plot(PP...; size=(900, 810), layout=(3, 3))
end
```

## ブートストラップ法

$n$ が奇数のとき, 一様分布 $\on{Uniform}(0, 1)$ のサイズ $n$ の標本の中央値の真の分布は $\Beta((n+1)/2, (n+1)/2)$ になる(下の方にある[二項分布とベータ分布の関係](#二項分布とベータ分布の関係)の節または[順序統計量 - Wikipedia](https://ja.wikipedia.org/wiki/%E9%A0%86%E5%BA%8F%E7%B5%B1%E8%A8%88%E9%87%8F) を参照せよ).

$n$ が10以上の偶数の場合には, 一様分布 $\on{Uniform}(0, 1)$ のサイズ $n$ の標本の中央値の真の分布(Dirichlet分布の積分で構成できる)の良い近似として, $n' = n+1$ のときの $\Beta((n'+1)/2, (n'+1)/2)$ を採用できる. $n$ が奇数の場合と違ってベータ分布中で使う $n$ を1増やしていることに注意せよ.  実際に良い近似になっていることについては以下のセルのグラフを見よ. グラフを見れば $n'=n+1$ を使った方が真の分布をよく近似しており, $n\ge 10$ で近似の精度が十分に高そうなことがわかる. 

```julia
plot_mediandist_approx()
```

$n$ が奇数のときは $n'=n$ とおき, $n$ が偶数のときには $n'=n+1$ とおいて, 標本サイズ $n$ に対して, ベータ分布 $\on{beta}$ を次のように定める:

$$
\on{beta}=\Beta((n'+1)/2, (n'+1)/2)
$$

このとき, 信頼係数 $1-\alpha$ の中央値の信頼区間 $[L, U]$ を次のように構成できる(ブートストラップ法):

$$
\begin{aligned}
&
L = \quantile(X, \quantile(\on{beta}, \alpha/2)),
\\ &
U = \quantile(X, \quantile(\on{beta}, 1-\alpha/2)).
\end{aligned}
$$

対応するP値函数は次のように書ける:

$$
\pval_{\on{bootstrap}}(X, a) = \min\left(
\begin{array}{l}
1 \\
2\cdf(\on{beta}, \ecdf(X)(a)) \\
2(1 - \cdf(\on{beta}, \ecdf(X)(a))) \\
\end{array}
\right).
$$

この方法で定義されたP値函数では第一種の過誤の確率が名目有意水準よりも大きくなってしまう場合が出て来ることに注意する必要がある.

信頼区間とP値の概念は表裏一体である. このことについては

* 竹内啓『数理統計学』 p.103
* 竹村彰通『現代数理統計学』 p.202
* 久保川達也『現代数理統計学の基礎』 p.169

などの教科書を参照せよ. P値函数と信頼区間の対応は, 与えられたデータについて, P値函数の値がα以上になるパラメータの範囲が信頼区間に一致するという条件で与えられる.

データが与えられたとき, 横軸をパラメータとするP値函数のグラフを描いたとき, そのグラフを高さ $\alpha$ で切断して得られる区間が信頼係数 $1-\alpha$ の信頼区間になる.

異なる方法で構成された信頼区間を比較するには, 対応するP値函数を比較すればよい.

```julia
"""(0,1)区間の一様分布のサイズnの標本の中央値の分布の近似ベータ分布(nが奇数ならばexact)."""
function beta_median(n)
    n += iseven(n)
    Beta((n+1)/2, (n+1)/2)
end

"""ブートストラップ法による中央値の信頼区間"""
function ci_median_bootstrap(X::AbstractVector; α = 0.05)
    beta = beta_median(length(X))
    L = quantile(X, quantile(beta, α/2))
    U = quantile(X, quantile(beta, 1 - α/2))
    L, U
end

"""ブートストラップ法で使う標本の経験分布の標本の中央値の累積分布函数"""
function cdf_median_bootstrap(X::AbstractVector, a)
    beta = beta_median(length(X))
    cdf(beta, ecdf(X)(a))
end

"""`pval_median_bootstrap(X, a)` は `X` が何であっても `cdf_median_bootstrap(X, a)` が実装されていれば自動的に使えるようになる."""
function pval_median_bootstrap(X, a)
    c = cdf_median_bootstrap(X, a)
    min(1, 2c, 2(1 - c))
end
```

```julia
Random.seed!(3734649)

# テストサンプルの生成
dist = Gamma(2, 3)
n = 40
X = rand(dist, n)

# 信頼区間の計算
L, U = ci_median_bootstrap(X; α = 0.05)
@show [L, U]

# プロット
Q1 = histogram(X; norm=true, alpha=0.3, bin=-1:2:21, label="data")
vline!([median(X)]; label="median of data", lw=1.5, c=2, ls=:dash)
plot!([L, U], zeros(2); label="bootstrap ci", lw=10, c=2)
plot!(dist, -1, 21; label="true dist", c=:blue)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")

Q2 = plot(x -> pval_median_bootstrap(X, x), -1, 21; label="P-value")
vline!([median(X)]; label="median of data", lw=1.5, c=2, ls=:dash)
plot!([L, U], fill(0.05, 2); label="bootstrap ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

plot(Q1, Q2; size=(800, 300))
```

## 二項分布を用いた確率の正確な計算に帰着する方法

連続的母集団分布の中央値(真の中央値)を $a$ と書くとき, その母集団分布のサイズ $n$ の標本中の値で $a$ 以下のものがちょうど $k$ 個になる確率は二項分布 $\Binomial(n, 1/2)$ で値 $k$ が生じる確率に等しい.  この事実を使えば二項検定の場合と同様にして, 信頼区間とP値函数を定義できる. 

標本 $X=(X_1,\ldots,X_n)$ を小さな順に並べたもの(sortしたもの)を $X(1)\le\cdots\le X(n)$ と書く.

二項分布 $\on{bin}$ を　$\on{bin} = \Binomial(n, 1/2)$ と定める.

このとき, 信頼係数 $1-\alpha$ の中央値の信頼区間 $[L, U]$ を次のように構成できる:

$$
\begin{aligned}
&
L = X(\quantile(\on{bin}, \alpha/2)),
\\ &
U = X(\quantile(\on{bin}, 1 - \alpha/2) + 1).
\end{aligned}
$$

ここで $\quantile(\on{bin}, p)$ は二項分布 $\on{bin}$ の累積分布函数 $\cdf(\on{bin}, j)$ の値が $p$ 以下になる最大の $j$ になり, $X(i)$ は上で定義したように標本中の値で下から $i$ 番目に小さなものである.

対応するP値函数は次のように書ける:

$$
\pval_{\on{binomial}}(X, a) = \min\left(
\begin{array}{l}
1 \\
2\cdf(\on{bin}, k) \\
2(1 - \cdf(\on{bin}, k-1)) \\
\end{array}
\right).
$$

ここで $k$ は標本中の値で $a$ 以下のものの個数である. すなわち $k$ は $X(i) \le a$ となる最大の $i$ に等しい.

この方法で定義されたP値函数では第一種の過誤の確率を確実に名目有意水準以下にできる. しかし, その分だけ検出力は落ちてしまう可能性がある.

信頼区間とP値の概念は表裏一体である. このことについては

* 竹内啓『数理統計学』 p.103
* 竹村彰通『現代数理統計学』 p.202
* 久保川達也『現代数理統計学の基礎』 p.169

などの教科書を参照せよ. P値函数と信頼区間の対応は, 与えられたデータについて, P値函数の値がα以上になるパラメータの範囲が信頼区間に一致するという条件で与えられる.

データが与えられたとき, 横軸をパラメータとするP値函数のグラフを描いたとき, そのグラフを高さ $\alpha$ で切断して得られる区間が信頼係数 $1-\alpha$ の信頼区間になる.

異なる方法で構成された信頼区間を比較するには, 対応するP値函数を比較すればよい.

```julia
"""サイズnの標本中の連続的母集団分布の中央値以下になる値の個数の分布"""
bin_median(n) = Binomial(n, 1/2)

"""二項分布に帰着して作った中央値の信頼区間"""
function ci_median_binomial(X::AbstractVector; α = 0.05)
    bin = bin_median(length(X))
    X′ = sort(X)
    L = X′[quantile(bin, α/2)]
    U = X′[quantile(bin, 1 - α/2) + 1]
    L, U
end

"""正確二項検定のP値に帰着して作った中央値のP値"""
function pval_median_binomial(X::AbstractVector, a)
    bin = bin_median(length(X))
    k  = count(≤(a), X)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end
```

```julia
Random.seed!(3734649)

# テストサンプルの生成
dist = Gamma(2, 3)
n = 40
X = rand(dist, n)

# 信頼区間の計算
@show ci_bst = ci_median_bootstrap(X; α = 0.05)
@show ci_bin = ci_median_binomial(X; α = 0.05)

# プロット
Q3 = histogram(X; norm=true, alpha=0.3, bin=-1:2:21, label="data")
vline!([median(X)]; label="median of data", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bin), zeros(2); label="binomial ci", lw=10, c=2)
plot!(dist, -1, 21; label="true dist", c=:blue)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")

Q4 = plot(x -> pval_median_binomial(X, x), -1, 21; label="binomial P-value")
vline!([median(X)]; label="median of data", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bin), fill(0.05, 2); label="binomial ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

plot(Q3, Q4; size=(800, 300))
```

```julia
@show ci_bst = ci_median_bootstrap(X; α = 0.05)
@show ci_bin = ci_median_binomial(X; α = 0.05)
println()
flush(stdout)

plot(Q2, Q4; size=(800, 300)) |> display
println()
flush(stdout)

x = range(2, 10, 401)
plot(x, x -> pval_median_bootstrap(X, x); label="bootstrap P-value")
plot!(x, x -> pval_median_binomial(X, x); label="binomial P-value")
plot!(collect(ci_bst), fill(0.056, 2); label="bootstrap ci", lw=4, c=1)
plot!(collect(ci_bin), fill(0.044, 2); label="binomial ci", lw=4, c=2)
plot!(; xtick=-1:21, ytick=[0:0.05:0.1; 0.2:0.1:1])
title!("$(name(dist)), n=$n")
```

```julia
Random.seed!(3734649)

# テストサンプルの生成
dist = Gamma(2, 3)
n = 41
X = rand(dist, n)

# 信頼区間の計算
@show ci_bst = ci_median_bootstrap(X; α = 0.05)
@show ci_bin = ci_median_binomial(X; α = 0.05)
println()
flush(stdout)

# プロット
R2 = plot(x -> pval_median_bootstrap(X, x), -1, 21; label="bootstrap P-value")
vline!([median(X)]; label="median of data", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bst), fill(0.05, 2); label="bootstrap ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

R4 = plot(x -> pval_median_binomial(X, x), -1, 21; label="binomial P-value")
vline!([median(X)]; label="median of data", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bin), fill(0.05, 2); label="binomial ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

plot(R2, R4; size=(800, 300)) |> display
println()
flush(stdout)

x = range(2, 10, 401)
plot(x, x -> pval_median_bootstrap(X, x); label="bootstrap P-value")
plot!(x, x -> pval_median_binomial(X, x); label="binomial P-value")
plot!(collect(ci_bst), fill(0.056, 2); label="bootstrap ci", lw=4, c=1)
plot!(collect(ci_bin), fill(0.044, 2); label="binomial ci", lw=4, c=2)
plot!(; xtick=-1:21, ytick=[0:0.05:0.1; 0.2:0.1:1])
title!("$(name(dist)), n=$n")
```

## 二項分布とベータ分布の関係

二項分布 $\Beta(n, p)$ における累積分布函数はベータ分布の累積分布函数で記述できる:

$$
\begin{aligned}
&
\sum_{j=0}^k \binom{n}{j} p^j(1-p)^{n-j} = \int_p^1 \frac{t^k(1-k)^{n-k-1}}{B(k+1, n-k)}\,dt,
\\ &
\sum_{j=k}^n \binom{n}{j} p^j(1-p)^{n-j} = \int_0^p \frac{t^{k-1}(1-k)^{n-k}}{B(k, n-k+1)}\,dt.
\end{aligned}
$$

以下のセルに数値的な確認がある. 証明は両辺を $p$ で微分しても得られるし, 右辺の積分で部分積分を繰り返しても得られる. もしくはより確率論的に, 一様分布 $\on{Uniform}(0,1)$ のサイズ $n$ の標本中の $p$ 以下の数値の個数が $k$ 以上になる確率を2通りに記述することによって, 後者の公式の直観的な説明を得ることもできる. 右辺のベータ分布における確率は, 一様分布 $\on{Uniform}(0,1)$ のサイズ $n$ の標本中の下から $k$ 番目に小さな値が $p$ 以下になる確率として得られる.  一様分布 $\on{Uniform}(0,1)$ のサイズ $n$ の標本中の下から $k$ 番目に小さな値が $p$ 以下になる確率の密度と $dt$ の積は, $n$ 個を $k-1$ 個, $1$ 個, $n-k$ 個に分割する方法の個数と $k-1$ 個が $t$ 以下になる確率と $1$ 個が $t$ を含む微小区間に含まれる確率 $dt$ と $n-k$ 個が $t$ より大きくなる確率の積

$$
\begin{aligned}
\frac{n!}{(k-1)!1!(n-k)!}\, t^{k-1}\,dt\,(1-t)^{n-k} &=
\frac{\Gamma(n+1)}{\Gamma(k)\Gamma(n-k+1)}\, t^{k-1}\,dt\,(1-t)^{n-k} \\ &=
\frac{1}{B(k, n-k+1)}\,t^{k-1}(1-k)^{n-k}\,dt
\end{aligned}
$$

になる. これを $t$ について $0$ から $p$ まで積分したものが, 上の後者の公式の右辺になる. 以上の説明は, 一様分布 $\on{Uniform}(0,1)$ のサイズ $n$ の標本の順序統計量(標本中の $k$ 番目に小さな値)の分布がベータ分布で表されることの説明にもなっている.

コンピュータの実装での対応する公式は以下のようになる:

$$
\begin{aligned}
&
\cdf(\Binomial(n, p), k) = \on{ccdf}(\Beta(k+1, n-k), p), 
\\ &
\on{ccdf}(\Binomial(n, p), k-1) = \cdf(\Beta(k, n-k+1), p).
\end{aligned}
$$

以上の関係を使うと, 二項分布による記述をベータ分布による記述に書き換えたり, パラメータが整数のベータ分布による記述を二項分布による記述に書き換えたりできる. 以上の中央値の区間推定や検定の記述はこの書き換えによって見かけ上異なる記述が得られる. 例えば, 

* 奥村晴彦, [ヒストグラムから中央値・分位数とその信頼区間を求める](https://oku.edu.mie-u.ac.jp/~okumura/python/hist-median.html)

に書いてある中央値の区間推定と検定の式とこのノートの式を比較するときには注意が必要である.

```julia
n, p = 10, 0.4
@show [cdf(Binomial(n, p), k) for k in 0:n-1] .== [ccdf(Beta(k+1, n-k), p) for k in 0:n-1]
@show [ccdf(Binomial(n, p), k-1) for k in 1:n] .== [cdf(Beta(k, n-k+1), p) for k in 1:n]

plot(; legend=:bottomright)
plot!(0:n-1, [cdf(Binomial(n, p), k) for k in 0:n-1]; label="cdf(Binomial(n, p), k)")
plot!(0:n-1, [ccdf(Beta(k+1, n-k), p) for k in 0:n-1]; label="ccdf(Beta(k+1, n-k), p)", ls=:dash)
```

## 第一種の過誤の確率の比較

```julia
function sim_pval_median(; dist = Gamma(2, 3), n = 40, L = 10^5)
    a = median(dist)
    pval_bst = Vector{Float64}(undef, L)
    pval_bin = Vector{Float64}(undef, L)
    tmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, tmp[Threads.threadid()])
        pval_bst[i] = pval_median_bootstrap(X, a)
        pval_bin[i] = pval_median_binomial(X, a)
    end
    pval_bst, pval_bin
end

function plot_probtype1error(; dist = Gamma(2, 3), n = 40, L = 10^5)
    pval_bst, pval_bin = sim_pval_median(; dist, n, L)
    ecdf_bst = ecdf(pval_bst)
    ecdf_bin = ecdf(pval_bin)
    
    α = range(0, 1, 401)
    P1 = plot(; legend=:bottomright)
    plot!(α, α -> ecdf_bst(α); label="bootstrap")
    plot!(α, α -> ecdf_bin(α); label="binomial", ls=:dash)
    plot!([0, 1], [0, 1]; label="", ls=:dot, c=:black)
    plot!(; xtick=0:0.1:1, ytick=0:0.1:1)
    title!("$(name(dist)), n=$n")
    plot!(; xlabel="nominal significance level α", ylabel = "probability of type I error")

    α = range(0, 0.1, 401)
    P2 = plot(; legend=:bottomright)
    plot!(α, α -> ecdf_bst(α); label="bootstrap")
    plot!(α, α -> ecdf_bin(α); label="binomial", ls=:dash)
    plot!([0, 0.1], [0, 0.1]; label="", ls=:dot, c=:black)
    plot!(; xtick=0:0.01:1, ytick=0:0.01:1)
    title!("$(name(dist)), n=$n")
    plot!(; xlabel="nominal significance level α", ylabel = "probability of type I error")
    
    plot(P1, P2; size=(800, 400), leftmargin=3Plots.mm, bottommargin=3Plots.mm)
end
```

第一種の過誤の確率は母集団分布によらない.

```julia
for dist in (Normal(2, 3), Gamma(2, 3), Exponential(), LogNormal())
    plot_probtype1error(; dist, n = 20) |> display
    println(); flush(stdout)
end
```

```julia
for dist in (Normal(2, 3), Gamma(2, 3), Exponential(), LogNormal())
    plot_probtype1error(; dist, n = 21) |> display
    println(); flush(stdout)
end
```

第一種の過誤の確率は名目有意水準に近い方がよい. この点に関して2つの方法の優劣は付け難い.

```julia
dist = Uniform()
for n in (10, 20, 40, 80, 160, 320, 640)
    plot_probtype1error(; dist, n) |> display
    println(); flush(stdout)
end
```

```julia
dist = Uniform()
for n in (10, 20, 40, 80, 160, 320, 640)
    n += 1
    plot_probtype1error(; dist, n) |> display
    println(); flush(stdout)
end
```

## ヒストグラムの中央値の信頼区間


### ヒストグラムデータ

```julia
Random.seed!(3734649)

dist, n = Gamma(2, 3), 40
X = rand(dist, n)
bin = floor(minimum(X)):ceil(maximum(X))
H = fit(Histogram, X, bin)
@show H
@show H.weights
plot(H; alpha=0.3, label="size-100 sample")
```

### ヒストグラムから作られる分布

ヒストグラムから各ビンごとに一様分布に従う確率分布を作ることができる.

```julia
function histogramdist(h::Histogram)
    e = h.edges[1]
    w = h.weights
    u = [Uniform(e[i], e[i+1]) for i in eachindex(e)[1:end-1]]
    p = w/sum(w)
    MixtureModel(u, p)
end
```

```julia
Hdist = histogramdist(H)
@show Hdist
plot(x -> pdf(Hdist, x), (extrema(Hdist) .+ (-1, 1))...; label="histogram dist")
```

ヒストグラムから作られた確率分布 $\on{Hdist}$ から平均や分散や中央値などを計算できる.

```julia
@show mean(Hdist)
@show var(Hdist)
@show median(Hdist)
@show quantile.(Hdist, (0.25, 0.50, 0.75));
```

しかし, 対称性に問題がある.

```julia
hX1 = fit(Histogram, [1, 2, 4, 5], 0.5:5.5)
@show hX1dist = histogramdist(hX1)
@show mhX1 = median(hX1dist)
plot(hX1; alpha=0.3, label="", title="hX1", size=(400, 250))
vline!([mhX1]; label="", lw=3)
```

```julia
hX2 = fit(Histogram, [1, 3, 4, 5], 0.5:5.5)
@show hX2dist = histogramdist(hX2)
@show qhX2 = quantile(hX2dist, 0.25)
plot(hX2; alpha=0.3, label="", title="hX1", size=(400, 250))
vline!([qhX2]; label="", lw=3)
```

そこで, `prevfloat(p)` と `nextfloat(p)` での `quantile` の平均を取ることにする.

```julia
myquantile(X, p) = (p ≤ 0 || 1 ≤ p) ? quantile(X, p) : (quantile(X, prevfloat(p)) + quantile(X, nextfloat(p)))/2
mymedian(X) = myquantile(X, 0.5)
myquantile(H::Histogram, p) = myquantile(histogramdist(H), p)
```

```julia
@show mymhX1 = mymedian(hX1dist)
plot(hX1; alpha=0.3, label="", title="hX1", size=(400, 250))
vline!([mymhX1]; label="", lw=3)
```

```julia
@show myqhX2 = myquantile(hX2dist, 0.25)
plot(hX2; alpha=0.3, label="", title="hX2", size=(400, 250))
vline!([myqhX2]; label="", lw=3)
```

### ブートストラップ法

```julia
function ci_median_bootstrap(H::Histogram; α = 0.05)
    Hdist = histogramdist(H)
    beta = beta_median(sum(H.weights))
    L = myquantile(Hdist, quantile(beta, α/2))
    U = myquantile(Hdist, quantile(beta, 1 - α/2))
    L, U
end

function cdf_median_bootstrap(H::Histogram, a)
    Hdist = histogramdist(H)
    n = sum(H.weights)
    beta = beta_median(n)
    cdf(beta, cdf(Hdist, a))
end
```

```julia
ci_median_bootstrap(H)
```

`pval_ian_bootstrap` メソッドは `cdf_median_bootstrap` メソッドが実装されれば自動的に使えるようになる.

```julia
pval_median_bootstrap(H, 3.8)
```

### 二項分布に帰着

$n$ が奇数の場合には $n$ を1増やしておかないとブートストラップ法とのずれが大きくなる.

```julia
function ci_median_binomial(H::Histogram; α = 0.05)
    Hdist = histogramdist(H)
    n = round(Int, sum(H.weights))
    bin = bin_median(n)
    L = myquantile(Hdist, (quantile(bin,   α/2) - 0.5)/n)
    U = myquantile(Hdist, (quantile(bin, 1-α/2) + 0.5)/n)
    L, U
end

function pval_median_binomial(H::Histogram, a)
    Hdist = histogramdist(H)
    n = round(Int, sum(H.weights))
    bin = bin_median(n)
    c = cdf(Hdist, a)
    min(1, 2cdf(bin, n*c + 0.5), 2ccdf(bin, n*c - 0.5))
end
```

```julia
sum(H.weights)
```

```julia
ci_median_binomial(H)
```

```julia
pval_median_binomial(H, 3.8)
```

### 比較

```julia
# 信頼区間の計算
@show ci_bst = ci_median_bootstrap(H; α = 0.05)
@show ci_bin = ci_median_binomial(H; α = 0.05)
println()
flush(stdout)

# プロット
P2 = plot(x -> pval_median_bootstrap(H, x), 2, 15; label="bootstrap P-value")
vline!([median(Hdist)]; label="median of histogram", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bst), fill(0.05, 2); label="bootstrap ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

P3 = plot(x -> pval_median_binomial(H, x), 2, 15; label="binomial P-value")
vline!([median(Hdist)]; label="median of histogram", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bin), fill(0.05, 2); label="binomial ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

plot(P2, P3; size=(800, 300)) |> display
println()
flush(stdout)

x = range(2, 9, 401)
plot(x, x -> pval_median_bootstrap(H, x); label="bootstrap P-value")
plot!(x, x -> pval_median_binomial(H, x); label="binomial P-value")
plot!(collect(ci_bst), fill(0.056, 2); label="bootstrap ci", lw=4, c=1)
plot!(collect(ci_bin), fill(0.044, 2); label="binomial ci", lw=4, c=2)
plot!(; xtick=-1:21, ytick=[0:0.05:0.1; 0.2:0.1:1])
title!("$(name(dist)), n=$n, bin=$bin")
```

```julia
Random.seed!(3734649)

dist, n = Gamma(2, 3), 41
X = rand(dist, n)
bin = floor(minimum(X)):ceil(maximum(X))
H = fit(Histogram, X, bin)
@show H
@show H.weights
plot(H; alpha=0.3, label="size-100 sample")
```

```julia
Hdist = histogramdist(H)

# 信頼区間の計算
@show ci_bst = ci_median_bootstrap(H; α = 0.05)
@show ci_bin = ci_median_binomial(H; α = 0.05)
println()
flush(stdout)

# プロット
P2 = plot(x -> pval_median_bootstrap(H, x), 2, 15; label="bootstrap P-value")
vline!([median(Hdist)]; label="median of histogram", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bst), fill(0.05, 2); label="bootstrap ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

P3 = plot(x -> pval_median_binomial(H, x), 2, 15; label="binomial P-value")
vline!([median(Hdist)]; label="median of histogram", lw=1.5, c=2, ls=:dash)
plot!(collect(ci_bin), fill(0.05, 2); label="binomial ci", lw=4, c=2)
vline!([median(dist)]; label="true median", lw=1.5, c=:blue, ls=:dashdot)
title!("$(name(dist)), n=$n")
plot!(; ytick=[0:0.05:0.1; 0.2:0.1:1])

plot(P2, P3; size=(800, 300)) |> display
println()
flush(stdout)

x = range(2, 9, 401)
plot(x, x -> pval_median_bootstrap(H, x); label="bootstrap P-value")
plot!(x, x -> pval_median_binomial(H, x); label="binomial P-value")
plot!(collect(ci_bst), fill(0.056, 2); label="bootstrap ci", lw=4, c=1)
plot!(collect(ci_bin), fill(0.044, 2); label="binomial ci", lw=4, c=2)
plot!(; xtick=-1:21, ytick=[0:0.05:0.1; 0.2:0.1:1])
title!("$(name(dist)), n=$n, bin=$bin")
```

### ヒストグラムを経由した場合の第一種の過誤の確率

```julia
function sim_pval_median_hist(; dist = Gamma(2, 3), n = 40, L = 10^5)
    a = median(dist)
    pval_bst = Vector{Float64}(undef, L)
    pval_bin = Vector{Float64}(undef, L)
    tmp = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        X = rand!(dist, tmp[Threads.threadid()])
        bin = floor(minimum(X)):ceil(maximum(X))
        H = fit(Histogram, X, bin)
        pval_bst[i] = pval_median_bootstrap(H, a)
        pval_bin[i] = pval_median_binomial(H, a)
    end
    pval_bst, pval_bin
end

function plot_probtype1error_hist(; dist = Gamma(2, 3), n = 40, L = 10^5)
    pval_bst, pval_bin = sim_pval_median_hist(; dist, n, L)
    ecdf_bst = ecdf(pval_bst)
    ecdf_bin = ecdf(pval_bin)
    
    α = range(0, 1, 401)
    P1 = plot(; legend=:bottomright)
    plot!(α, α -> ecdf_bst(α); label="bootstrap")
    plot!(α, α -> ecdf_bin(α); label="binomial", ls=:dash)
    plot!([0, 1], [0, 1]; label="", ls=:dot, c=:black)
    plot!(; xtick=0:0.1:1, ytick=0:0.1:1)
    title!("$(name(dist)), n=$n")
    plot!(; xlabel="nominal significance level α", ylabel = "probability of type I error")

    α = range(0, 0.1, 401)
    P2 = plot(; legend=:bottomright)
    plot!(α, α -> ecdf_bst(α); label="bootstrap")
    plot!(α, α -> ecdf_bin(α); label="binomial", ls=:dash)
    plot!([0, 0.1], [0, 0.1]; label="", ls=:dot, c=:black)
    plot!(; xtick=0:0.01:1, ytick=0:0.01:1)
    title!("$(name(dist)), n=$n")
    plot!(; xlabel="nominal significance level α", ylabel = "probability of type I error")
    
    plot(P1, P2; size=(800, 400), leftmargin=3Plots.mm, bottommargin=3Plots.mm)
end

function plot_probtype1error_hist_iter(; dist = Uniform(0, √(12 * 18)))
    for n in (10, 20, 40, 80, 160, 320, 640)
        plot_probtype1error_hist(; dist, n) |> display
        println(); flush(stdout)
    end
    for n in (10, 20, 40, 80, 160, 320, 640)
        n += 1
        plot_probtype1error_hist(; dist, n) |> display
        println(); flush(stdout)
    end
end
```

```julia
uniform = Uniform(0, round(√(12 * 18); digits=3))
normal = Normal(2, √18)
gamma = Gamma(2, 3)
exponential = Exponential(round(√18; digits=3))
lognormal = LogNormal(0, round(√log((1 + √(1+4*18))/2); digits=3))

@show var(uniform)
@show var(normal)
@show var(gamma)
@show var(exponential)
@show var(lognormal);
```

```julia
for dist in (uniform, normal, gamma, exponential, lognormal)
    plot_probtype1error_hist(; dist, n = 21) |> display
    println(); flush(stdout)
end
```

```julia
plot_probtype1error_hist_iter(; dist = uniform)
```

```julia
plot_probtype1error_hist_iter(; dist = normal)
```

```julia
plot_probtype1error_hist_iter(; dist = gamma)
```

```julia
plot_probtype1error_hist_iter(; dist = exponential)
```

```julia
plot_probtype1error_hist_iter(; dist = lognormal)
```

## SARS-CoV-2の変異株B.1.1.529系統（オミクロン株）の潜伏期間の推定：暫定報告

* [SARS-CoV-2の変異株B.1.1.529系統（オミクロン株）の潜伏期間の推定：暫定報告](https://www.niid.go.jp/niid/ja/2019-ncov/2551-cepr/10903-b11529-period.html)

より:

>データ２では、アルファ株症例1118例、オミクロン株症例113例が解析の対象となった。アルファ株症例の潜伏期間の中央値は3.4日（95％信頼区間：3.3-3.6）、オミクロン株症例は2.9日（95％信頼区間：2.5-3.2）であった。感染曝露から95％、99％が発症するまでの日数は、アルファ株症例ではそれぞれ8.7日、11.9日、オミクロン株症例ではそれぞれ7.1日、9.7日であった。

<img src="omi_per_f2.png" width="500">!

>感染曝露からの経過日数ごとの累積発症確率を表１に示す。アルファ株では10日目までに97.35％が発症するのに対して、オミクロン株では99.18%が発症すると推定された。

<img src="omi_per_t1.png" width="500">


### ヒストグラムから症例数を目で読み取る.

以下はグラフから読み取った症例数.

```julia
data = [
     1 184 19
     2 177 28
     3 231 31
     4 183 16
     5 131 11
     6  64  5
     7  71  2
     8  19  1
     9  32  0
    10  13  0
    11   6  0
    12   4  0
    13   1  0
    14   2  0
]

df = DataFrame(data, ["days after exposure", "Alpha n", "Omicron n"])
```

### ヒストグラムデータの中央値の信頼区間などを計算


「$n$ 日目に発症」を「暴露されてから $n-0.5$ 日目の始めから $n+0.5$ 日目の直前までのあいだに発症」と解釈してヒストグラムに変換し, 中央値と中央値の信頼区間を計算してみる.

```julia
bin = 0.5:14.5
@show hist_Alpha = Histogram((bin,), df[!, "Alpha n"], :right, false)
@show hist_Omicron = Histogram((bin,), df[!, "Omicron n"], :right, false);
```

```julia
@show size_Alpha = sum(df[!, "Alpha n"])
@show median_Alpha = mymedian(hist_Alpha)
@show ci_median_bootstrap(hist_Alpha; α=0.05)
@show ci_median_binomial(hist_Alpha; α=0.05)
@show ci_median_bootstrap(hist_Alpha; α=0.01)
@show ci_median_binomial(hist_Alpha; α=0.01)
println()
@show size_Omicron = sum(df[!, "Omicron n"])
@show median_Alpha = mymedian(hist_Omicron)
@show ci_median_bootstrap(hist_Omicron; α=0.05)
@show ci_median_binomial(hist_Omicron; α=0.05)
@show ci_median_bootstrap(hist_Omicron; α=0.01)
@show ci_median_binomial(hist_Omicron; α=0.01);

println("\n"*"-"^78*"\n")

@show ci_Alpha = ci_median_binomial(hist_Alpha; α=0.05)
@show ci_Omicron = ci_median_binomial(hist_Omicron; α=0.05)
println()

P1 = plot(hist_Alpha; alpha=0.3, label="")
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
plot!(collect(ci_Alpha), fill(-0.01*maximum(hist_Alpha.weights), 2); lw=5, c=:blue, label="95% CI of median")
vline!([mymedian(hist_Alpha)]; c=:blue, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(hist_Alpha, 0.95)]; c=:blue, label="", ls=:dash)
vline!([myquantile(hist_Alpha, 0.99)]; c=:blue, label="", ls=:dash)
title!("Alpha n = $size_Alpha")

P2 = plot(hist_Omicron; alpha=0.3, label="", c=2)
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
plot!(collect(ci_Omicron), fill(-0.01*maximum(hist_Omicron.weights), 2); lw=5, c=:red, label="95% CI of median")
vline!([mymedian(hist_Omicron)]; c=:red, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(hist_Omicron, 0.95)]; c=:red, label="", ls=:dash)
vline!([myquantile(hist_Omicron, 0.99)]; c=:red, label="", ls=:dash)
title!("Omicron n = $size_Omicron")

plot(P1, P2; size=(800, 220))
```

中央値の信頼区間だけを計算しても面白くないのだが, このノート内では一般のquantileの信頼区間は実装していないのでプロットしなかった. 原理的には一般のquantileの信頼区間の実装も易しい.


[SARS-CoV-2の変異株B.1.1.529系統（オミクロン株）の潜伏期間の推定：暫定報告](https://www.niid.go.jp/niid/ja/2019-ncov/2551-cepr/10903-b11529-period.html)にあった以下のグラフに近いものが得られた. ただし, 以上のプロットはガンマ分布モデルでフィッティングせずに, ヒストグラムのデータから直接計算した結果である.

<img src="omi_per_f2.png" width="500">


### ガンマ分布, 対数正規分布, Weibull分布でフィッティング

ガンマ分布, 対数正規分布, Weibull分布などでフィッティングしたバージョンも作ってみよう.

```julia
sample_Alpha = foldl(vcat, (fill(i, df[i, "Alpha n"]) for i in 1:14))
sample_Omicron = foldl(vcat, (fill(i, df[i, "Omicron n"]) for i in 1:14))

@show gamma_Alpha = fit_mle(Gamma, sample_Alpha)
@show gamma_Omicron = fit_mle(Gamma, sample_Omicron)

P1 = plot(hist_Alpha; alpha=0.3, label="")
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
vline!([mymedian(gamma_Alpha)]; c=:blue, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(gamma_Alpha, 0.95)]; c=:blue, label="", ls=:dash)
vline!([myquantile(gamma_Alpha, 0.99)]; c=:blue, label="", ls=:dash)
title!("Alpha n = $size_Alpha / Gamma fitting")
plot!(x -> size_Alpha * pdf(gamma_Alpha, x); label="", c=:blue, ls=:dot, lw=1.5)

P2 = plot(hist_Omicron; alpha=0.3, label="", c=2)
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
vline!([mymedian(gamma_Omicron)]; c=:red, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(gamma_Omicron, 0.95)]; c=:red, label="", ls=:dash)
vline!([myquantile(gamma_Omicron, 0.99)]; c=:red, label="", ls=:dash)
title!("Omicron n = $size_Omicron / Gamma fitting")
plot!(x -> size_Omicron * pdf(gamma_Omicron, x); label="", c=:red, ls=:dot, lw=1.5)

plot(P1, P2; size=(800, 220))
```

```julia
"""https://www.niid.go.jp/niid/ja/2019-ncov/2551-cepr/10903-b11529-period.html"""
data_gamma_orig = [
1 6.29 8.55
2 23.1 30.41
3 42.42 53.05
4 59.46 70.69
5 72.67 82.65
6 82.16 90.12
7 88.63 94.53
8 92.90 97.04
9 95.63 98.43
10 97.35 99.18
11 98.41 99.57
12 99.05 99.78
13 99.44 99.89
14 99.67 99.94
]

df[!, "Alpha cdf orig"] = data_gamma_orig[:, 2]
df[!, "Alpha cdf"] = round.(100cdf.(gamma_Alpha, 1:14); digits=2)
df[!, "Omicron cdf orig"] = data_gamma_orig[:, 3]
df[!, "Omicron cdf"] = round.(100cdf.(gamma_Omicron, 1:14); digits=2)
df
```

```julia
sample_Alpha = foldl(vcat, (fill(i, df[i, "Alpha n"]) for i in 1:14))
sample_Omicron = foldl(vcat, (fill(i, df[i, "Omicron n"]) for i in 1:14))

lognormal_Alpha = fit_mle(LogNormal, sample_Alpha)
lognormal_Omicron = fit_mle(LogNormal, sample_Omicron)

size_Alpha = sum(df[!, "Alpha n"])
size_Omicron = sum(df[!, "Omicron n"])

P1 = plot(hist_Alpha; alpha=0.3, label="")
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
vline!([mymedian(lognormal_Alpha)]; c=:blue, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(lognormal_Alpha, 0.95)]; c=:blue, label="", ls=:dash)
vline!([myquantile(lognormal_Alpha, 0.99)]; c=:blue, label="", ls=:dash)
title!("Alpha n = $size_Alpha / LogNormal fitting")
plot!(x -> size_Alpha * pdf(lognormal_Alpha, x); label="", c=:blue, ls=:dot, lw=1.5)

P2 = plot(hist_Omicron; alpha=0.3, label="", c=2)
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
vline!([mymedian(lognormal_Omicron)]; c=:red, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(lognormal_Omicron, 0.95)]; c=:red, label="", ls=:dash)
vline!([myquantile(lognormal_Omicron, 0.99)]; c=:red, label="", ls=:dash)
title!("Omicron n = $size_Omicron / LogNormal fitting")
plot!(x -> size_Omicron * pdf(lognormal_Omicron, x); label="", c=:red, ls=:dot, lw=1.5)

plot(P1, P2; size=(800, 220))
```

```julia
sample_Alpha = foldl(vcat, (fill(i, df[i, "Alpha n"]) for i in 1:14))
sample_Omicron = foldl(vcat, (fill(i, df[i, "Omicron n"]) for i in 1:14))

weibull_Alpha = fit_mle(Weibull, sample_Alpha)
weibull_Omicron = fit_mle(Weibull, sample_Omicron)

size_Alpha = sum(df[!, "Alpha n"])
size_Omicron = sum(df[!, "Omicron n"])

P1 = plot(hist_Alpha; alpha=0.3, label="")
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
vline!([mymedian(weibull_Alpha)]; c=:blue, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(weibull_Alpha, 0.95)]; c=:blue, label="", ls=:dash)
vline!([myquantile(weibull_Alpha, 0.99)]; c=:blue, label="", ls=:dash)
title!("Alpha n = $size_Alpha / Weibull fitting")
plot!(x -> size_Alpha * pdf(weibull_Alpha, x); label="", c=:blue, ls=:dot, lw=1.5)

P2 = plot(hist_Omicron; alpha=0.3, label="", c=2)
plot!(; xlim=(-0.5, 14.5), xtick=0:14)
vline!([mymedian(weibull_Omicron)]; c=:red, label="median, 95%, 99%", ls=:dash)
vline!([myquantile(weibull_Omicron, 0.95)]; c=:red, label="", ls=:dash)
vline!([myquantile(weibull_Omicron, 0.99)]; c=:red, label="", ls=:dash)
title!("Omicron n = $size_Omicron / Weibull fitting")
plot!(x -> size_Omicron * pdf(weibull_Omicron, x); label="", c=:red, ls=:dot, lw=1.5)

plot(P1, P2; size=(800, 220))
```

```julia
@show aic_gamma_Alpha = -2loglikelihood(gamma_Alpha, sample_Alpha) + 4
@show aic_lognormal_Alpha = -2loglikelihood(lognormal_Alpha, sample_Alpha) + 4
@show aic_weibull_Alpha = -2loglikelihood(weibull_Alpha, sample_Alpha) + 4
println()
@show aic_gamma_Omicron = -2loglikelihood(gamma_Omicron, sample_Omicron) + 4
@show aic_lognormal_Omicron = -2loglikelihood(lognormal_Omicron, sample_Omicron) + 4
@show aic_weibull_Omicron = -2loglikelihood(weibull_Omicron, sample_Omicron) + 4;
```

確かに, ガンマ分布モデル, 対数正規分布モデル, Weibull分布モデルの3つの中ではガンマ分布モデルのAICの値が一番小さい.

```julia

```
