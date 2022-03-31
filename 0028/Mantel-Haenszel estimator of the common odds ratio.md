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

# 共通オッズ比のMantel-Haenszelの推定量と最尤推定量の関係

* 黒木玄
* 2022-03-29～2022-03-31
$
\newcommand\op{\operatorname}
$

## 文献紹介

共通オッズ比のMantel-Haenszelの推定量と最尤推定量の関係については次の文献に書いてある:

* Yanagimoto, Yakemi and Eiji Yamamoto, Eiji. Simple linear approximations to the likelihood equation for combining evidence in multiple 2×2 tables: A critique of conventional procedures. Ann. Inst. Statist. Math., Vol. 37, 1985, 37-49. \[[link](https://link.springer.com/article/10.1007/BF02481079)\] \[[pdf](https://www.ism.ac.jp/editsec/aism/pdf/037_1_0037.pdf)\]

この文献ではMantel-Haenszel以外の共通オッズ比の推定量も扱っている.  共通オッズ比のMantel-Haenszelの推定量の対数の分散のRobins-Breslow-Greenlandの推定量の解説が次の文献に書いてある:

* Silcocks, Paul.  An easy approach to the Robins-Breslow-Greenland variance estimator. Epidemiol Perspect Innov. 2005. \[[link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1270683/https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1270683/)\] \[[pdf](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1270683/pdf/1742-5573-2-9.pdf)\]


## 最尤法とMantel-Haenszelの共通オッズ比の推定量の関係


### 統計モデル

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

の正規分布に近似的に従うことを示せる. (二項分布に関する中心極限定理と所謂デルタ法(逆数を取る操作の一次近似)の簡単な応用で示せる.)

以下では $\omega_1 = \cdots = \omega_K = \omega$ が成立していると仮定する.  $\omega$ を __共通オッズ比__ と呼ぶ.

この仮定のもとで以上で扱っている統計モデルの独立なパラメータは $q_1,\ldots,q_K, \omega$ の $K+1$ 個になる.


### 共通オッズ比のMantel-Haenszelの推定量

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

大数の法則より, $m_k, n_k$ 達が大きなとき, 近似

$$
\hat\omega_k =
\frac{(a_k/m_k)(d_k/n_k)}{(b_k/m_k)(c_k/n_k)}\approx
\frac{p_k(1-q_k)}{(1-p_k)q_k} =
\omega_k = 
\omega
$$

が成立しているので, 近似

$$
\hat\omega_{\op{MH}} \approx
\frac
{\sum_{k=1}^K \hat{w}_k \omega}
{\sum_{k=1}^K \hat{w}_k} = \omega.
$$

が成立する.  すなわち大数の法則による近似が有効な場合にはMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ は共通オッズ比 $\omega$ を近似している.

これはMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ が共通オッズ比 $\omega$ の一致推定量になっていることを意味している.

共通オッズ比が $1$ の特殊な場合 ($\omega = 1$ の場合)にはさらに良いことを言える.

大数の法則より, 近似

$$
\frac{1}{\hat{w}_k} = \frac{m_k + n_k}{b_k c_k} =
\frac{1}{(b_k/m_k)(c_k/n_k)}\left(\frac{1}{m_k} + \frac{1}{n_k}\right) \approx
\frac{1}{p_k (1 - q_k)}\left(\frac{1}{m_k} + \frac{1}{n_k}\right)
$$

が成立しているので,  $\omega_k = \omega = 1$ (すなわち $p_k = q_k$) ならば,

$$
\frac{1}{\hat{w}_k} \approx
\frac{1}{p_k (1 - p_k)}\left(\frac{1}{m_k} + \frac{1}{n_k}\right) = 
\omega^2\left(\frac{1}{m_k p_k (1 - p_k)} + \frac{1}{n_k q_k (1 - q_k)}\right) = v_k.
$$

すなわち $\hat{w}_k$ の逆数は $\hat\omega_k = (a_k d_k)/(b_k c_k)$ の漸近的な分散 $v_k$ を近似する.

一般に, 正の実数達 $\sigma_1^2,\ldots,\sigma_K^2 > 0$ が与えられていて,  統計モデルが共通平均 $\mu$ と既知の分散 $\sigma_k^2$ を持つ $K$ 個の正規分布の積であるとき, 共通平均 $\mu$ の共通平均 $\mu$ の標本 $X_1,\ldots,X_n$ に対する最尤推定量が

$$
\hat\mu =
\frac
{\sum_{k=1}^K X_k/\sigma_k^2}
{\sum_{k=1}^K 1/\sigma_k^2}
$$

になることを簡単な計算で示せる. これは分散の逆数を重みとする荷重平均になっている.  分散が小さな $X_k$ の方が確率的に共通平均の真の値に近くなる傾向があるはずなので, そのような $X_k$ を大きな重みで足し上げた方が共通平均の良い推定量になりそうである.  この荷重平均は実際にそのようになっている.  (この荷重平均は $\sigma_1^2 = \cdots = \sigma_K^2 = \sigma^2$ の場合の標本平均 $\bar{X}=\frac{1}{K}\sum_{k=1}^K X_k$ の一般化になっている.) 

上で示した共通オッズ比が $1$ の場合 ($\omega_k = \omega = 1$) におけるMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ に関する結果は, 共通オッズ比が $1$ の特殊な場合には, Mantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ が $\mu=\omega=1$, $\sigma_k^2=v_k$ の場合の 漸近的最尤推定量 $\hat\mu$ を近似していることを意味している.

ただし, これは共通のオッズ比が $1$ の特殊な場合におけるMantel-Haenszelの推定量 $\hat\omega_{\op{MH}}$ と最尤推定量の関係に過ぎない.


### 共通オッズ比の最尤推定量

$K$ 個の独立な2×2の分割表の各々が2つの二項分布の積分布に従っているという統計モデルのパラメータ空間を共通オッズ比 $\omega$ を持つ場合に制限して得られる統計モデルを考えている. 独立なパラメータは $q_1,\ldots,q_K,\omega$ の $K+1$ 個になる.  (ロジスティックモデルで書き直すこともできる.)

このモデルにおける標本

$$
A_k = \begin{bmatrix}
a_k & b_k \\
c_k & d_k \\
\end{bmatrix}
\quad (k = 1,\ldots,K)
$$

の統計モデルのパラメータ達の最尤推定量 $\hat{q}_1, \ldots, \hat{q}_K, \hat\omega$ を以下のようにして計算できることを示せる. (対数尤度函数を偏微分したものが $0$ になるという方程式を整理すればよい.)

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

これを共通オッズ比の __最尤方程式__ と呼ぶことにする. さらに

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


![IMG_4660.PNG](attachment:30191889-742c-468d-bb7d-179086c2e979.PNG)


### 最尤方程式の線形近似とMantel-Haenszel推定量の関係

この節の内容は上で引用した Yamagimoto-Yamamoto (1985) の Proposition 1 の紹介である.

$y = \delta_k(x)$ は $x = \hat\omega_k = (a_k d_k)/(b_k c_k)$ で $y = 0$ となり, $x=1$ で $y = (a_k d_k - b_k c_k)/N_k$ となる.  その2点を通る直線

$$
y =
\frac{(a_k d_k - b_k c_k)/N_k}{1 - (a_k d_k)/(b_k c_k)}\left(x - \frac{a_k d_k}{b_k c_k}\right) =
-\frac{b_k c_k}{N_k}x + \frac{a_k d_k}{N_k}
$$

で $y = \delta_k(x)$ を近似することにしよう(これはかなり大胆な近似である).

この近似を採用すると, 最尤方程式 ($*$) を次の一次方程式 ($**$) で近似することになる:

$$
-\left(\sum_{k=1}^K \frac{b_k c_k}{N_k}\right)x + \sum_{k=1}^K \frac{a_k d_k}{N_k} = 0.
\tag{$**$}
$$

これを共通オッズ比の __線形近似最尤方程式__ と呼ぶことにする. この線形近似最尤方程式の解はMantel-Haenszel推定量

$$
\hat\omega_{\op{MH}} =
\frac
{\sum_{k=1}^K a_k d_k/N_k}
{\sum_{k=1}^K b_k c_k/N_k}.
$$

に一致する.  このように共通オッズ比のMantel-Haenszel推定量は最尤方程式($*$)を線形近似最尤方程式($**$)で近似することによって得られる.


![IMG_4659.PNG](attachment:6526d398-c43a-4063-bcc9-6561e3b84311.PNG)


### Mantel-Haenszelの推定量の対数の分散のRobins-Breslow-Greenlandの推定量

この節の内容は上で引用した Silcocks (2005) の紹介になっている.

Mantel-Haenszelの推定量は次のように表される:

$$
\hat\omega_{\op{MH}} =
\frac{\sum_{k=1}^K a_k d_k/N_k}{\sum_{k=1}^K b_k c_k/N_k} =
\frac{\sum_{k=1}^K \hat{w}_k \hat\omega_k}{\sum_{k=1}^K \hat{w}_k}, \quad
\hat\omega_k = \frac{a_k d_k}{b_k c_k}, \quad
\hat{w}_k = \frac{b_k c_k}{N_k}
$$

デルタ法(対数の一次近似)によって,

$$
\op{var}(\log \hat\omega_{\op{MH}}) \approx
\frac{\op{var}(\hat\omega_{\op{MH}})}{\hat\omega_{\op{MH}}^2}.
$$


$k$ 番目の分割表のオッズ比の最尤推定量 $\hat\omega_k$ は近似的に平均が共通オッズ比 $\omega$ で分散が

$$
\op{var}(\hat\omega_k) \approx
\omega^2\left(
\frac{1}{m_k p_k} +
\frac{1}{m_k (1 - p_k)} +
\frac{1}{n_k q_k} +
\frac{1}{n_k (1 - q_k)}
\right)
$$

の正規分布に近似的に従うことを示せるのであった. さらに, 

$$
\op{var}(\hat\omega_k) \approx \hat\omega_{\op{MH}}^2\left(\frac{1}{a_k}+\frac{1}{b_k}+\frac{1}{c_k}+\frac{1}{d_k}\right)
$$

という近似も使えると仮定し, 重み達 $\hat{w}_k$ が定数であるかのような近似が成立していると仮定する.  このとき,

$$
\op{var}(\hat\omega_{\op{MH}}) \approx
\frac{\sum_{k=1}^K \hat{w}_k^2 \op{var}(\hat\omega_k)}{\left(\sum_{k=1}^K \hat{w}_k\right)^2}.
$$

以上の近似をまとめると次の近似が得られる:

$$
\op{var}(\log\hat\omega_{\op{MH}}) \approx
\frac{\op{var}(\hat\omega_{\op{MH}})}{\hat\omega_{\op{MH}}^2} \approx
\frac
{\sum_{k=1}^K \hat{w}_k^2\left(\frac{1}{a_k}+\frac{1}{b_k}+\frac{1}{c_k}+\frac{1}{d_k}\right)}
{\left(\sum_{k=1}^K \hat{w}_k\right)^2}.
$$


さらに, $\hat\omega_{\op{MH}} \approx (a_k d_k)/(b_k c_k)$ という近似も使えると仮定する. このとき

$$
\frac{1}{a_k}+\frac{1}{b_k}+\frac{1}{c_k}+\frac{1}{d_k} =
\frac{a_k + d_k}{a_k d_k}+\frac{b_k + c_k}{b_k c_k} \approx
\frac{1}{b_k c_k}\left(\frac{a_k + d_k}{\hat\omega_{\op{MH}}}+ b_k + c_k\right).
$$

これと $\hat{w}_k = b_k c_k/N_k$ を使うと, 上で示した近似式より,

$$
\op{var}(\log\hat\omega_{\op{MH}}) \approx
\frac
{\sum_{k=1}^K b_k c_k (\hat\omega_{\op{MH}}^{-1}(a_k + d_k)+ b_k + c_k) / N_k^2}
{\left(\sum_{k=1}^K b_k c_k/N_k\right)^2} =
\frac
{\sum_{k=1}^K \hat\omega_{\op{MH}}^{-1} b_k c_k (a_k + d_k + \hat\omega_{\op{MH}}(b_k + c_k)) / N_k^2}
{\left(\sum_{k=1}^K b_k c_k/N_k\right)^2}.
$$

$(a_k, b_k)$ と $(c_k, d_k)$ を交換すると, $\hat\omega_{\op{MH}}$ は逆数になり,  $\log\hat\omega_{\op{MH}}$ は $-1$ 倍になって分散は変わらないので,

$$
\op{var}(\log\hat\omega_{\op{MH}}) \approx
\frac
{\sum_{k=1}^K a_k d_k (\hat\omega_{\op{MH}}(b_k + c_k) + a_k + d_k) / N_k^2}
{\left(\sum_{k=1}^K a_k d_k/N_k \right)^2}.
$$


記号の簡単のため, ここで次のようにおく:

$$
P_k = \frac{a_k + d_k}{N_k}, \quad
Q_k = \frac{b_k + c_k}{N_k}, \quad
R_k = \frac{a_k d_k}{N_k}, \quad
S_k = \frac{b_k c_k}{N_k}, \quad
R = \sum_{k=1}^K R_k, \quad
S = \sum_{k=1}^K S_k.
$$

このとき, 上で示した2つの近似式の平均を取ると,

$$
\begin{aligned}
\op{var}(\log\hat\omega_{\op{MH}}) & \approx
\frac
{\sum_{k=1}^K (R^2\hat\omega_{\op{MH}}^{-1} b_k c_k + S^2 a_k d_k)
(a_k + d_k + \hat\omega_{\op{MH}}(b_k + c_k)) / N_k^2}
{2 R^2 S^2}
\\ & =
\frac
{\sum_{k=1}^K (R^2\hat\omega_{\op{MH}}^{-1} S_k + S^2 R_k)
(P_k + \hat\omega_{\op{MH}}Q_k)}
{2 R^2 S^2}
\end{aligned}
$$


さらに分子分母を $S^2$ で割って, $\hat\omega_{\op{MH}} = R/S$ を使うと, 

$$
\begin{aligned}
\op{var}(\log\hat\omega_{\op{MH}}) & \approx
\frac
{\sum_{k=1}^K (\hat\omega_{\op{MH}} S_k + R_k) (P_k + \hat\omega_{\op{MH}}Q_k)}
{2 R^2}
\\ & =
\frac
{\sum_{k=1}^K (P_k R_k + \hat\omega_{\op{MH}}(P_k S_k + Q_k R_k) + \hat\omega_{\op{MH}}^2 Q_k S_k)}
{2 R^2}
\\ & =
\frac{\sum_{k=1}^K P_k R_k}{2R^2} +
\frac{\sum_{k=1}^K (P_k S_k + Q_k R_k)}{2RS} +
\frac{\sum_{k=1}^K Q_k S_k}{2S^2}.
\end{aligned}
$$

この最後の量を __Mantel-Haenszelの推定量の対数の分散のRobins-Breslow-Greenlandの推定量__ と呼び, $V_{\op{RBG}}$ と書くことにする. 

Mantel-Haenszelの推定量の対数 $\log\hat\omega_{\op{MH}}$ は平均 $\log \omega$ (共通オッズ比の対数), 分散 $V_{\op{RBG}}$ の正規分布に近似的に従う.

この事実を使って共通オッズ比に関するP値函数や信頼区間を構成することができる.

以下の数値例を見ればわかるように, このRBG公式は非常に良い公式であり, 最尤法(スコア検定)によって構成したP値函数をRBG公式を使って構成したP値函数は近似的によく一致する.


## 数値例

```julia
using RCall
using Distributions
using StatsPlots
default(fmt = :png)
using Roots
```

```julia
safediv(x, y) = x == 0 ? x/one(y) : x/y
safesqrt(x) = √max(0, x)

oddsratio(a, b, c, d) = safediv(a*d, b*c)
oddsratio(A::AbstractVecOrMat) = oddsratio(A...)

function delta(a, b, c, d, ω)
    A, B, C = 1 - ω, a + d + ω*(b + c), a*d - ω*b*c
    safediv(2C, B + safesqrt(B^2 - 4A*C))
end
delta(A, ω) = delta(A..., ω)

function linear_approx_delta(a, b, c, d, ω)
    N = a + b + c + d
    -(b*c/N)*ω + a*d/N
end
linear_approx_delta(A::AbstractVecOrMat, ω) = linear_approx_delta(A..., ω)

eachmatrix(A::AbstractArray{<:Any, 3}) = eachslice(A; dims=3)

function maximum_likelihood_equation(A::AbstractArray{<:Any, 3}, x)
    sum(A -> delta(A, x), eachmatrix(A))
end

function linear_approx_maximum_likelihood_equation(A::AbstractArray{<:Any, 3}, x)
    sum(A -> linear_approx_delta(A, x), eachmatrix(A))
end

function mantel_haenszel_estimator(A::AbstractArray{<:Any, 3})
    num = sum(A -> A[1,1]*A[2,2]/sum(A), eachmatrix(A))
    den = sum(A -> A[1,2]*A[2,1]/sum(A), eachmatrix(A))
    num/den
end

function maximum_linkelihood_estimator(A::AbstractArray{<:Any, 3})
    f(t) = maximum_likelihood_equation(A, exp(t))
    logomegahat = find_zero(f, 0.0, Order2())
    omegahat = exp(logomegahat)
end
```

以下の2つのセルの内容については

* https://github.com/genkuroki/public/blob/main/0028/Mantel-Haenszel.ipynb

も参照せよ.

```julia
function vardelta(a, b, c, d, ω)
    N = a + b + c + d
    δ = delta(a, b, c, d, ω)
    (N - 1)/N/(1/(a - δ) + 1/(b + δ) + 1/(c + δ) + 1/(d - δ))
end
vardelta(A::AbstractVecOrMat, ω) = vardelta(A..., ω)

function chisq_mantel_haenszel(A::AbstractArray{<:Any, 3}, ω = 1.0)
    num = sum(A -> delta(A, ω), eachmatrix(A))^2
    den = sum(A -> vardelta(A, ω), eachmatrix(A))
    num/den
end

function pvalue_chisq(A::AbstractArray{<:Any, 3}, ω = 1.0)
    chisq = chisq_mantel_haenszel(A, ω)
    ccdf(Chisq(1), chisq)
end

function ci_chisq(A::AbstractArray{<:Any, 3}, α = 0.05)
    f(t) = pvalue_chisq(A, exp(t)) - α
    logci = find_zeros(f, -1e1, 1e1)
    exp(first(logci)), exp(last(logci))
end

function mh_chisq(A::AbstractArray{<:Any, 3}; ω₀ = 1.0, α = 0.05)
    common_odds_ratio = maximum_linkelihood_estimator(A)
    chisq = chisq_mantel_haenszel(A, ω₀)
    df = 1
    p_value = ccdf(Chisq(df), chisq)
    conf_int = ci_chisq(A, α)
    (; common_odds_ratio, ω₀, p_value, α, conf_int, chisq, df)
end
```

```julia
function mantel_haenszel_robins_breslow_greenland(A::AbstractArray{<:Any, 3})
    @views a, b, c, d  = A[1,1,:], A[1,2,:], A[2,1,:], A[2,2,:]
    abcd = zip(a, b, c, d)
    R = sum(((a, b, c, d),) -> a*d/(a+b+c+d), abcd)
    S = sum(((a, b, c, d),) -> b*c/(a+b+c+d), abcd)
    logOR = log(R) - log(S)
    sumPR = sum(((a, b, c, d),) -> (a+d)*a*d/(a+b+c+d)^2, abcd)
    sumPSplusQR = sum(((a, b, c, d),) -> ((a+d)*b*c + (b+c)*a*d)/(a+b+c+d)^2, abcd)
    sumQS = sum(((a, b, c, d),) -> (b+c)*b*c/(a+b+c+d)^2, abcd)
    SE² = sumPR/(2R^2) + sumPSplusQR/(2(R*S)) + sumQS/(2S^2)
    SE = √SE²
    (; logOR, SE)
end

function pvalue_mhrbg(A::AbstractArray{<:Any, 3}, ω = 1.0)
    (; logOR, SE) = mantel_haenszel_robins_breslow_greenland(A)
    normal = Normal(logOR, SE)
    min(1, 2cdf(normal, log(ω)), 2ccdf(normal, log(ω)))
end

function ci_mhrbg(A::AbstractArray{<:Any, 3}, α = 0.05)
    (; logOR, SE) = mantel_haenszel_robins_breslow_greenland(A)
    normal = Normal(logOR, SE)
    exp.(quantile.(normal, (α/2, 1 - α/2)))
end

function mh_rbg(A::AbstractArray{<:Any, 3}; ω₀ = 1.0, α = 0.05)
    (; logOR, SE) = mantel_haenszel_robins_breslow_greenland(A)
    common_odds_ratio = exp(logOR)
    normal = Normal(logOR, SE)
    p_value = min(1, 2cdf(normal, log(ω₀)), 2ccdf(normal, log(ω₀)))
    conf_int = exp.(quantile.(normal, (α/2, 1 - α/2)))
    (; common_odds_ratio, ω₀, p_value, α, conf_int)
end
```

### 例1

* https://twitter.com/mph_for_doctors/status/1226714150853849090

```julia
R"""
library(tidyverse)

df.forMH <- tribble(
  ~X, ~L, ~Y, ~n,
  0, 0, 0, 325,
  0, 0, 1, 273,
  0, 1, 0, 324,
  0, 1, 1, 363,
  1, 0, 0, 292,
  1, 0, 1, 321,
  1, 1, 0, 278,
  1, 1, 1, 323
)

df.MH1 <- uncount(df.forMH, weights = n)
"""

R"""
library(samplesizeCMH)

partial_tables <- table(df.MH1[,c("X","Y","L")])
partial_tables
"""
```

```julia
R"""
result_mh <- mantelhaen.test(partial_tables)
"""
```

```julia
R"""
mantelhaen.test(partial_tables)

library(jtools)
log.model <- glm(Y ~ X + L, data=df.MH1, family = binomial)
result_logistic <- summ(log.model, exp=T, digits=5)
"""
```

```julia
R"df_MH1 <- df.MH1"
@rget df_MH1
@rget partial_tables
```

```julia
A = partial_tables
@show A
@show mh_chisq(A)
@show mh_rbg(A)

f(x) = maximum_likelihood_equation(A, x)
g(x) = linear_approx_maximum_likelihood_equation(A, x)
plot(f, 0.5, 2; label="maximum likelihood equation")
plot!(g; label="Mantel-Haenszel linear approx.", ls=:dash)
hline!([0]; label="", c=:black, lw=0.5, ls=:dot)
plot!(; xlabel="common odds ratio", ylabel="score")
plot!(; xtick=0.1:0.1:3)
```

このグラフのMantel-Haenszel linear approx.が線形近似最尤方程式であり, その零点が共通オッズ比のMantel-Haenszel推定量になる.  最尤方程式(maximum likelihood equation)の零点が(unconditionalな場合の)最尤法で求めた共通オッズ比の推定量になる.

```julia
plot(title="P-value functions")
plot!(ω -> pvalue_chisq(A, ω), 0.7, 1.7; label="Mantel-Haenszel χ²")
plot!(ω -> pvalue_mhrbg(A, ω); label="Robins-Breslow-Greenland", ls=:dash)
plot!(; xtick=0.1:0.1:3, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
```

この場合に, 最尤法の場合のスコア検定のP値函数(上のグラフのMantel-Haenszel χ²)と共通オッズ比のMantel-Haenszel推定量の対数の分散のRobins-Breslow-Greenlandの推定量を使った正規分布近似で作ったP値函数(上のグラフのRobins-Breslow-Greenland)はほぼぴったり一致している.


### 例2

* https://twitter.com/MinatoNakazawa/status/1202358323409850369
    * http://minato.sip21c.org/im3r/20191204.html
        * http://minato.sip21c.org/epispecial/codes-for-Chapter8.R
            * Ten Studies: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(95)92163-X/fulltext Figure 2B
            * Eleventh Study: https://www.nejm.org/doi/full/10.1056/NEJM199810083391504 Table 3

```julia
TenStudies = [
    215 229 311-215 306-229
     38  33  59-38   51-33
    161 174 293-161 293-174
     76  88 164-76  163-88
    103 105 129-103 133-105
     65  67 120-65  125-67
     81  75 113-81  110-75
     48  63 160-48  159-63
     22  21  60-22  62-21
     56  51 137-56  140-51
]

ElevenStudies = [
    TenStudies
    468 480 697-468 685-480
]
```

```julia tags=[]
oddsratio.(eachrow(ElevenStudies))
```

```julia
A = A10 = reshape(TenStudies', 2, 2, :) |> collect

@show A
@show mh_chisq(A)
@show mh_rbg(A)

f(x) = maximum_likelihood_equation(A, x)
g(x) = linear_approx_maximum_likelihood_equation(A, x)
plot(f, 0.5, 1.5; label="maximum likelihood equation")
plot!(g; label="Mantel-Haenszel linear approx.", ls=:dash)
hline!([0]; label="", c=:black, lw=0.5, ls=:dot)
plot!(; xlabel="common odds ratio", ylabel="score")
plot!(; xtick=0.1:0.1:3)
```

```julia
plot(title="P-value functions of Ten Studies")
plot!(ω -> pvalue_chisq(A10, ω), 0.7, 1.1; label="Mantel-Haenszel χ²")
plot!(ω -> pvalue_mhrbg(A10, ω); label="Robins-Breslow-Greenland", ls=:dash)
plot!(; xtick=0.1:0.05:3, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
plot!(; xlabel="common odds ratio", ylabel="P-value")
```

```julia
A = A11 = reshape(ElevenStudies', 2, 2, :) |> collect

@show A
@show mh_chisq(A)
@show mh_rbg(A)

f(x) = maximum_likelihood_equation(A, x)
g(x) = linear_approx_maximum_likelihood_equation(A, x)
plot(f, 0.5, 1.5; label="maximum likelihood equation")
plot!(g; label="Mantel-Haenszel linear approx.", ls=:dash)
hline!([0]; label="", c=:black, lw=0.5, ls=:dot)
plot!(; xlabel="common odds ratio", ylabel="score")
plot!(; xtick=0.1:0.1:3)
```

```julia
plot(title="P-value functions")
plot!(ω -> pvalue_chisq(A11, ω), 0.7, 1.1; label="Mantel-Haenszel χ²")
plot!(ω -> pvalue_mhrbg(A11, ω); label="Robins-Breslow-Greenland", ls=:dash)
plot!(; xtick=0.1:0.05:3, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
plot!(; xlabel="common odds ratio", ylabel="P-value")
```

```julia
plot(title="P-value functions")
plot!(ω -> pvalue_chisq(A10, ω), 0.7, 1.1; label="Ten Studies")
plot!(ω -> pvalue_chisq(A11, ω); label="Eleven Studies", ls=:dash)
plot!(; xtick=0.1:0.05:3, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
plot!(; xlabel="common odds ratio", ylabel="P-value")
```

* https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(00)02163-2/fulltext

```julia
TenStudies = [
    274 276 311-274 306-276
     45  39  59-45   51-39
    215 222 293-215 293-222
    134 139 164-134 163-139
    119 127 129-119 133-127
     86  92 120-86  125-92
     81  75 113-81  110-75
     95  92 168-95  162-92
     51  54  60-54   62-54
     56  51 137-56  140-51
]

ElevenStudies = [
    TenStudies
    468 480 698-468 687-480
]
```

```julia tags=[]
oddsratio.(eachrow(ElevenStudies))
```

```julia
A = reshape(ElevenStudies', 2, 2, :) |> collect

@show A
@show mh_chisq(A)
@show mh_rbg(A)

f(x) = maximum_likelihood_equation(A, x)
g(x) = linear_approx_maximum_likelihood_equation(A, x)
plot(f, 0.5, 1.5; label="maximum likelihood equation")
plot!(g; label="Mantel-Haenszel linear approx.", ls=:dash)
hline!([0]; label="", c=:black, lw=0.5, ls=:dot)
plot!(; xlabel="common odds ratio", ylabel="score")
plot!(; xtick=0.1:0.1:3)
```

```julia
plot(title="P-value functions")
plot!(ω -> pvalue_chisq(A, ω), 0.7, 1.2; label="Mantel-Haenszel χ²")
plot!(ω -> pvalue_mhrbg(A, ω); label="Robins-Breslow-Greenland", ls=:dash)
plot!(; xtick=0.1:0.05:3, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
plot!(; xlabel="common odds ratio", ylabel="P-value")
```

### 例3

```julia
data = [
    10 20 50 100
    150 50 50 20
    50 50 50 50
    200 160 150 180
]
```

```julia
oddsratio.(eachrow(data))
```

```julia
A = reshape(data', 2, 2, :) |> collect

@show A
@show mh_chisq(A)
@show mh_rbg(A)

f(x) = maximum_likelihood_equation(A, x)
g(x) = linear_approx_maximum_likelihood_equation(A, x)
plot(f, 0.5, 2.2; label="maximum likelihood equation")
plot!(g; label="Mantel-Haenszel linear approx.", ls=:dash)
hline!([0]; label="", c=:black, lw=0.5, ls=:dot)
plot!(; xlabel="common odds ratio", ylabel="score")
plot!(; xtick=0.1:0.1:3)
```

```julia
plot(title="P-value functions")
plot!(ω -> pvalue_chisq(A, ω), 0.7, 2; label="Mantel-Haenszel χ²")
plot!(ω -> pvalue_mhrbg(A, ω); label="Robins-Breslow-Greenland", ls=:dash)
plot!(; xtick=0.1:0.1:3, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
plot!(; xlabel="common odds ratio", ylabel="P-value")
```

### 例4

```julia
data = [
    10 5 5 10
    10 5 10 20
    10 10 5 10
    10 20 5 10
]
```

```julia
oddsratio.(eachrow(data))
```

```julia
A = reshape(data', 2, 2, :) |> collect

@show A
@show mh_chisq(A)
@show mh_rbg(A)

f(x) = maximum_likelihood_equation(A, x)
g(x) = linear_approx_maximum_likelihood_equation(A, x)
plot(f, 0.8, 3.6; label="maximum likelihood equation")
plot!(g; label="Mantel-Haenszel linear approx.", ls=:dash)
hline!([0]; label="", c=:black, lw=0.5, ls=:dot)
plot!(; xlabel="common odds ratio", ylabel="score")
plot!(; xtick=0.2:0.2:5)
```

このグラフのMantel-Haenszel linear approx.が線形近似最尤方程式であり, その零点が共通オッズ比のMantel-Haenszel推定量になる.  最尤方程式(maximum likelihood equation)の零点が(unconditionalな場合の)最尤法で求めた共通オッズ比の推定量になる. この場合にそれらは互いに少しずれている.

```julia
plot(title="P-value functions")
plot!(ω -> pvalue_chisq(A, ω), 0.5, 7; label="Mantel-Haenszel χ²")
plot!(ω -> pvalue_mhrbg(A, ω); label="Robins-Breslow-Greenland", ls=:dash)
plot!(; xtick=0.5:0.5:10, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
plot!(; xlabel="common odds ratio", ylabel="P-value")
```

この場合には, 最尤法の場合のスコア検定のP値函数(上のグラフのMantel-Haenszel χ²)と共通オッズ比のMantel-Haenszel推定量の対数の分散のRobins-Breslow-Greenlandの推定量を使った正規分布近似で作ったP値函数(上のグラフのRobins-Breslow-Greenland)は互いに少しずれている.


### 例5

* Rothman-Greenland-Lash, Modern Epideomology, 3rd edition, 2008, p.331, Table 15-5 

```julia
data = [
    3 9 104 1059
    1 3   5   86
]
```

```julia
oddsratio.(eachrow(data))
```

```julia
A = reshape(data', 2, 2, :) |> collect

@show A
@show mh_chisq(A)
@show mh_rbg(A)

f(x) = maximum_likelihood_equation(A, x)
g(x) = linear_approx_maximum_likelihood_equation(A, x)
plot(f, 0.8, 5; label="maximum likelihood equation")
plot!(g; label="Mantel-Haenszel linear approx.", ls=:dash)
hline!([0]; label="", c=:black, lw=0.5, ls=:dot)
plot!(; xlabel="common odds ratio", ylabel="score")
plot!(; xtick=0.2:0.2:10)
```

```julia
plot(title="P-value functions")
plot!(ω -> pvalue_chisq(A, ω), 0.5, 16; label="Mantel-Haenszel χ²")
plot!(ω -> pvalue_mhrbg(A, ω); label="Robins-Breslow-Greenland", ls=:dash)
plot!(; xtick=0:30, ytick=0:0.05:1)
vline!([1]; label="", c=:black, lw=0.5)
plot!(; xlabel="common odds ratio", ylabel="P-value")
```

```julia

```
