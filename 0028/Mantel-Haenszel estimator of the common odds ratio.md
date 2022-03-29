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

# 共通オッズ比に関するMantel-Haenszelの推定量と最尤推定量の関係

* 黒木玄
* 2022-03-29
$
\newcommand\op{\operatorname}
$

## 文献紹介

共通オッズ比に関するMantel-Haenszelの推定量と最尤推定量の関係については次の文献に書いてある:

* Yanagimoto, Yakemi and Eiji Yamamoto, Eiji. Simple linear approximations to the likelihood equation for combining evidence in multiple 2×2 tables: A critique of conventional procedures. Ann. Inst. Statist. Math., Vol. 37, 1985, 37-49. \[[link](https://link.springer.com/article/10.1007/BF02481079https://link.springer.com/article/10.1007/BF02481079)\] \[[pdf](https://www.ism.ac.jp/editsec/aism/pdf/037_1_0037.pdf)\]

この文献ではMantel-Haenszel以外の共通オッズ比の推定量も扱っている.  共通オッズ比に関するMantel-Haenszelの推定量の対数の分散に関するRobins-Breslow-Greenlandの推定量の解説が次の文献に書いてある:

* Silcocks, Paul.  An easy approach to the Robins-Breslow-Greenland variance estimator. Epidemiol Perspect Innov. 2005. \[[link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1270683/https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1270683/)\] \[[pdf](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1270683/pdf/1742-5573-2-9.pdfhttps://www.ncbi.nlm.nih.gov/pmc/articles/PMC1270683/pdf/1742-5573-2-9.pdf)\]


## 統計モデル

簡単のため $K$ 個の独立な2×2の分割表の各々は2つの二項分布の積分布に従っているという統計モデルを考える.  $k$ 番目の分割表を行列

$$
A_k = \begin{bmatrix}
a_k & b_k \\
c_k & d_k \\
\end{bmatrix}
$$

と書き, 

$$
N_k = a_k + b_k + c_k + d_k
$$

とおく.  $a_k + b_k = m_k$ と $c_k + d_k = n_k$ は定数であり, $a_k$, $c_k$ はそれぞれ独立な二項分布 $\op{Binomial}(m_k, p_k)$, $\op{Binomial}(n_k, q_k)$ に従っているとする.  $k$ 番目のオッズ比 $\omega_k$ を

$$
\omega_k = \frac{p_k(1 - q_k)}{(1 - p_k)q_k}
$$

と定める.  $0<p_k<1$ かつ $0<q_k<1$ のとき $\omega_k = 1$ と $p_k = q_k$ は同値である.

以下では $K$ が固定されており, $m_k$, $n_k$ 達が十分に大きい場合を扱う.

このとき, $K$ 番目のオッズ比の推定量

$$
\hat\omega_k = \frac{a_k d_k}{b_k c_k}
$$

は平均が $\omega_k$ で分散が

$$
\begin{aligned}
v_k &= \omega^2\left(
\frac{1}{m_k p_k} +
\frac{1}{m_k (1 - p_k)} +
\frac{1}{n_k q_k} +
\frac{1}{n_k (1 - q_k)}
\right)
\\ &= \omega^2 \left(
\frac{1}{m_k p_k (1 - p_k)} +
\frac{1}{n_k q_k (1 - q_k)}
\right)
\end{aligned}
$$

の正規分布に近似的に従うことを示せる. (二項分布に関する中心極限定理と所謂デルタ法=1次までのTaylor展開の簡単な応用で示せる.)

以下では $\omega_1 = \cdots = \omega_K = \omega$ が成立していると仮定する.  $\omega$ を __共通オッズ比__ と呼ぶ.

この仮定のもとで以上で扱っている統計モデルの独立なパラメータは $q_1,\ldots,q_K, \omega$ の $K+1$ 個になる.


## 共通オッズ比に関するMantel-Haenszelの推定量

共通オッズ比の __Mantel-Haenszelの推定量__ $\hat\omega_{\op{MH}}$ を次のように定める:

$$
\hat\omega_{\op{MH}} =
\frac
{\sum_{k=1}^K a_k d_k/N_k}
{\sum_{k=1}^K b_k c_k/N_k}.
$$

これは次のように書き直される:

$$
\hat\omega_{\op{MH}} =
\frac
{\sum_{k=1}^K \hat{w}_k \hat\omega_k}
{\sum_{k=1}^K \hat{w}_k}.
$$

ここで

$$
\hat\omega_k = \frac{a_k d_k}{b_k c_k}, \quad
\hat{w}_k = \frac{b_k c_k}{N_k}
$$

大数の法則より, $m_k, n_k$ 達が十分大きなとき, 近似

$$
\hat\omega_k \approx \omega_k = \omega, \quad
$$

が成立しているので, 近似

$$
\hat\omega_{\op{MH}} \approx
\frac
{\sum_{k=1}^K \hat{w}_k \omega}
{\sum_{k=1}^K \hat{w}_k} = \omega.
$$

が成立している.  すなわち Mantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ は共通オッズ比 $\omega$ を近似している.

これはMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ が共通オッズ比 $\omega$ の一致推定量になっていることを意味している.

共通オッズ比が $1$ の特殊な場合 ($\omega = 1$ の場合)にはさらに良いことを言える.

大数の法則より, 

$$
\frac{1}{\hat{w}_k} = \frac{m_k + n_k}{b_k c_k} =
\frac{1}{(m_k b_k)(n_k c_k)}\left(\frac{1}{m_k} + \frac{1}{n_k}\right) \approx
\frac{1}{p_k (1 - q_k)}\left(\frac{1}{m_k} + \frac{1}{n_k}\right)
$$

なので $\omega_k = \omega = 1$ (すなわち $p_k = q_k$) ならば,

$$
\begin{aligned}
\frac{1}{\hat{w}_k} &\approx
\frac{1}{p_k (1 - p_k)}\left(\frac{1}{m_k} + \frac{1}{n_k}\right) \\ &= 
\omega^2\left(\frac{1}{m_k p_k (1 - p_k)} + \frac{1}{n_k q_k (1 - q_k)}\right) = v_k.
\end{aligned}
$$

すなわち $\hat{w}_k$ の逆数は $\hat\omega_k = (a_k d_k)/(b_k c_k)$ の漸近的な分散 $v_k$ を近似する.

一般に, 与えられた $\sigma_1^2,\ldots,\sigma_K^2 > 0$ について,  統計モデルが共通の平均 $\mu$ と分散 $\sigma_k^2$ を持つ $K$ 個の正規分布の積であるとき, 標本 $X_1,\ldots,X_n$ に関する共通の平均 $\mu$ の最尤推定量が

$$
\hat\mu =
\frac
{\sum_{k=1}^K X_k/\sigma_k^2}
{\sum_{k=1}^K 1/\sigma_k^2}
$$

になることを簡単な計算で示せる. これは分散の逆数を重みとする荷重平均になっている.  分散が小さな $X_k$ の方が共通平均の真の値に近くなる傾向があるので, そのような $X_k$ を大きな重みで足し上げた方が良い推定量になりそうである.  この荷重平均は実際にそのようになっている.  (この荷重平均は $\sigma_1^2 = \cdots = \sigma_K^2 = \sigma^2$ の場合の標本平均 $\bar{X}=\frac{1}{K}\sum_{k=1}^K X_k$ の一般化になっている.) 

上で示した共通オッズ比が $1$ の場合 ($\omega_k = \omega = 1$) に関するMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ に関する結果は, 共通オッズ比が $1$ の場合にはMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ が近似的に漸近的な最尤推定量 $\hat\mu$ の形をしていることを意味している.

ただし, これは共通のオッズ比が $1$ の特殊な場合におけるMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ と最尤推定量の関係に過ぎない.


## 最尤推定量

$K$ 個の独立な2×2の分割表の各々が2つの二項分布の積分布に従っているという統計モデルのパラメータ空間を共通オッズ比 $\omega$ を持つ場合に制限して得られる統計モデルを考えている. 独立なパラメータは $q_1,\ldots,q_K,\omega$ の $K+1$ 個になる.  (ロジスティックモデルで書き直すこともできる.)

このモデルにおける標本

$$
A_k = \begin{bmatrix}
a_k & b_k \\
c_k & d_k \\
\end{bmatrix}
\quad (k = 1,\ldots,K)
$$

の最尤推定量 $\hat{q}_1, \ldots, \hat{q}_K, \hat\omega$ を以下のようにして計算できることを示せる. (対数尤度函数を偏微分したものが $0$ になるという方程式を整理すればよい.)

(1) 次の条件を満たす $\delta_k = \delta_k(\omega)$ を求める:

$$
\frac{(a_k - \delta_k)(d_k - \delta_k)}{(b_k + \delta_k)(c_k + \delta_k)} = \omega, \quad
-\min(b_k, c_k) < \delta_k < \min(a_k, d_k)
$$

具体的に $\delta_k = \delta_k(\omega)$ は2次方程式を解いて次のように表される:

$$
\delta_k = \delta_k(\omega) =
\frac{2C_k}{B_k + \sqrt{B_k^2 - 4A_k C_k}}.
$$

ここで

$$
A_k = 1 - \omega, \quad
B_k = a_k + d_k + \omega(b_k + c_k), \quad
C_k = a_k d_k - \omega b_k c_k.
$$

このとき, $\delta_k = \delta_k(\omega)$ は $\omega$ の単調減少函数になり, 次を満たしている:

$$
\delta_k(0) = \min(a_k, d_k), \quad
\delta_k(1) = \frac{a_k d_k - b_k c_k}{N_k}, \quad
\delta_k(\hat\omega_k) = 0, \quad
\delta_k(\infty) = -\min(b_k, c_k).
$$

$\hat\omega_k = (a_k d_k)/(b_k c_k)$ と定めたのであった.

(2) $x > 0$ に関する次の方程式の解を $\hat\omega$ とする:

$$
\sum_{k=1}^K \delta_k(x) = 0.
\tag{$*$}
$$

さらに

$$
\begin{bmatrix}
\hat{a}_k & \hat{b}_k \\
\hat{c}_k & \hat{d}_k \\
\end{bmatrix} =
\begin{bmatrix}
a_k - \delta_k(\hat\omega) & b_k + \delta_k(\hat\omega) \\
c_k + \delta_k(\hat\omega) & d_k - \delta_k(\hat\omega) \\
\end{bmatrix}
$$

とおき, $\hat{q}_k = \hat{c}_k/n_k$ と定める.  このとき, $p_k$ の最尤推定量は $\hat{p}_k = \hat{a}_k/m_k$ になる.

以上における方程式($*$)を __最尤方程式__ と呼ぶことにする.


## 最尤方程式の線形近似とMantel-Haenszel推定量

この節の内容は上で引用した Yamagimoto-Yamamoto (1985) の Proposition 1 の紹介である.

$y = \delta_k(x)$ は $x = \hat\omega_k = (a_k d_k)/(b_k c_k)$ で $y = 0$ となり, $x=1$ で $y = (a_k d_k - b_k c_k)/N_k$ となる.  その2点を通る直線

$$
y =
\frac{(a_k d_k - b_k c_k)/N_k}{1 - (a_k d_k)/(b_k c_k)}\left(x - \frac{a_k d_k}{b_k c_k}\right) =
-\frac{b_k c_k}{N_k}x + \frac{a_k d_k}{N_k}
$$

で $y = \delta_k(x)$ を近似することにしよう(これはかなり大胆な近似である).  このとき最尤方程式 ($*$) を次の1次方程式 ($**$) で近似することになる:

$$
-\left(\sum_{k=1}^K \frac{b_k c_k}{N_k}\right)x + \sum_{k=1}^K \frac{a_k d_k}{N_k} = 0.
\tag{$**$}
$$

これを __線形近似最尤方程式__ と呼ぶことにする. この線形近似最尤方程式の解はMantel-Haenszel推定量

$$
\hat\omega_{\op{MH}} =
\frac
{\sum_{k=1}^K a_k d_k/N_k}
{\sum_{k=1}^K b_k c_k/N_k}.
$$

に一致する.  このようにMantel-Haenszel推定量は最尤方程式($*$)を線形近似最尤方程式($**$)で近似することによって得られる.

```julia
using Distributions
using StatsPlots
default(fmt = :png)
```

```julia
oddsratio(a, b, c, d) = a*d/(b*c)
oddsratio(A) = oddsratio(A...)

function delta(a, b, c, d, ω)
    ω == 1 && return (a*d - b*c)/(a+b+c+d)
    A, B, C = 1 - ω, a + d + ω*(b + c), a*d - ω*b*c
    2C/(B + √(B^2 - 4A*C))
end

delta(A, ω) = delta(A..., ω)

function linapprox(a, b, c, d, ω)
    N = a + b + c + d
    -(b*c/N)*ω + a*d/N
end

linapprox(A, ω) = linapprox(A..., ω)
```

```julia
A = [
     5 10
    10  5
]

ω = range(0, 2, 1000)
plot(ω, ω -> delta(A, ω))
plot!(ω, ω -> linapprox(A, ω))
```

```julia
A = [
    100   50
     50  100
]

ω = range(0, 6, 1000)
plot(ω, ω -> delta(A, ω))
plot!(ω, ω -> linapprox(A, ω))
```

```julia
delta(A, 1), delta(A, oddsratio(A))
```

```julia
linapprox(A, 1), linapprox(A, oddsratio(A))
```

```julia

```
