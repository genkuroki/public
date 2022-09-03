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
    display_name: Julia 1.8.0
    language: julia
    name: julia-1.8
---

# 正規分布モデルの共役事前分布によるベイズ統計

* 黒木玄
* 2022-09-03
$
\newcommand\ds{\displaystyle}
\newcommand\op[1]{{\operatorname{#1}}}
\newcommand\R{{\mathbb R}}
\newcommand\var{\op{var}}
\newcommand\cov{\op{cov}}
\newcommand\ybar{{\bar y}}
\newcommand\sigmahat{{\hat\sigma}}
$

<!-- #region toc=true -->
<h1>目次<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#正規分布モデルの共役事前分布とその応用" data-toc-modified-id="正規分布モデルの共役事前分布とその応用-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>正規分布モデルの共役事前分布とその応用</a></span><ul class="toc-item"><li><span><a href="#逆ガンマ正規分布" data-toc-modified-id="逆ガンマ正規分布-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>逆ガンマ正規分布</a></span></li><li><span><a href="#Bayes更新" data-toc-modified-id="Bayes更新-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Bayes更新</a></span></li><li><span><a href="#μの周辺事前分布" data-toc-modified-id="μの周辺事前分布-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>μの周辺事前分布</a></span></li><li><span><a href="#事前・事後予測分布" data-toc-modified-id="事前・事後予測分布-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>事前・事後予測分布</a></span></li><li><span><a href="#Jeffreys事前分布の場合" data-toc-modified-id="Jeffreys事前分布の場合-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>Jeffreys事前分布の場合</a></span></li><li><span><a href="#通常の信頼区間と予測区間との比較" data-toc-modified-id="通常の信頼区間と予測区間との比較-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>通常の信頼区間と予測区間との比較</a></span></li><li><span><a href="#データの数値から事前分布を決めた場合" data-toc-modified-id="データの数値から事前分布を決めた場合-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>データの数値から事前分布を決めた場合</a></span></li></ul></li></ul></div>
<!-- #endregion -->

```julia
using Distributions
using StatsPlots
default(fmt=:png, size=(400, 250),
    titlefontsize=10, tickfontsize=6, guidefontsize=9,
    plot_titlefontsize=10)
using SymPy
```

## 正規分布モデルの共役事前分布とその応用


### 逆ガンマ正規分布

平均 $\mu\in\R$, 分散 $v=\sigma>0$ の正規分布の確率密度函数を次のように表す:

$$
p_\op{N}(y|\mu, v) =
\frac{1}{\sqrt{2\pi v}}\exp\left(-\frac{1}{2v}(y-\mu)^2\right)
\quad (y\in \R).
$$

パラメータ $\kappa, \theta > 0$ の逆ガンマ分布の確率密度函数を次のように表す:

$$
p_\op{IG}(v|\kappa,\theta) =
\frac{\theta^\kappa}{\Gamma(\kappa)}
v^{-\kappa-1}\exp\left(-\frac{\theta}{v}\right)
\quad (v > 0).
$$

$v$ がこの逆ガンマ分布に従う確率変数だとすると, 

$$
\begin{aligned}
&
\frac{1}{v} \sim
\op{Gamma}\left(\kappa,\, \frac{1}{\theta}\right) =
\frac{1}{2\theta}\op{Gamma}\left(\frac{2\kappa}{2},\, 2\right) =
\frac{1}{2\theta}\op{Chisq}(2\kappa),
\\ &
E[v] = \frac{\theta}{\kappa - 1}, \quad
\var(v) = \frac{E[v]^2}{\kappa - 2}.
\end{aligned}
$$

逆ガンマ正規分布の密度函数を次のように定義する:

$$
\begin{aligned}
p_\op{IGN}(\mu,v|\mu_*, v_*, \kappa, \theta) &=
p_\op{N}(\mu|\mu_*, v_* v) p_\op{IG}(v|\kappa, \theta)
\\ &\propto
v^{-(\kappa-1/2)-1}
\exp\left(-\frac{1}{v}\left(\theta + \frac{1}{2v_*}(\mu-\mu_*)^2\right)\right).
\end{aligned}
$$

この逆ガンマ正規分布の密度函数に従う確率変数を $\mu,v$ と書くと,

$$
E[v] = \frac{\theta}{\kappa-1}, \quad
\var(v) = \frac{E[v]^2}{\kappa-2}, \quad
\cov(\mu, v) = 0, \quad
E[\mu] = \mu_*, \quad
\var(\mu) = v_* E[v].
$$

逆ガンマ分布は正規分布の共役事前分布になっている.


### Bayes更新

以下において, $A$ と $B$ が $\mu, v$ に関する定数倍を除いて等しいことを $A\propto B$ と書く.

データの数値 $y_1,\ldots,y_n$ が与えられたとき, 正規分布モデルの尤度函数は

$$
\prod_{i=1}^n p_\op{N}(y_i|\mu,v) \propto
v^{-n/2}\exp\left(-\frac{1}{2v}\sum_{i=1}^n(y_i-\mu)^2\right)
$$

の形になる. このとき,

$$
\ybar = \frac{1}{n}\sum_{i=1}^n y_i, \quad
\sigmahat^2 = \frac{1}{n}\sum_{i=1}^n (y_i - \ybar)^2.
$$

とおくと, 

$$
\sum_{i=1}^n(y_i-\mu)^2 = n(\mu - \ybar)^2 + n\sigmahat^2
$$

なので, 尤度を最大化する $\mu, v$ は $\mu=\ybar$, $v=\sigmahat^2$ になることがわかる.

さらに, 次が成立することもわかる:

$$
\begin{aligned}
&
\prod_{i=1}^n p_\op{N}(y_i|\mu,v)\times p_\op{IGN}(\mu,v|\mu_*, v_*, \kappa, \theta)
\\ &\propto
v^{-n/2}\exp\left(-\frac{1}{2v}(y_i-\mu)^2\right)\times
v^{-(\kappa+1/2)-1}
\exp\left(-\frac{1}{v}\left(\theta + \frac{1}{2v_*}(\mu-\mu_*)^2\right)\right)
\\ &\propto
v^{-(\kappa+n/2+1/2)-1}
\exp\left(-\frac{1}{v}\left(
\theta + \frac{n\sigmahat^2}{2} +
\frac{1+nv_*}{2v_*}\left(\mu - \frac{\mu_*+nv_*\ybar}{1+nv_*}\right)
\right)\right).
\end{aligned}
$$

ゆえに共役事前分布から得られる事後分布のパラメータは次のようになる:

$$
\begin{alignedat}{2}
&
\tilde\kappa = \kappa + \frac{n}{2} =
\frac{n}{2}\left(1 + \frac{2\kappa}{n}\right), \qquad
& &
\tilde\theta = \theta + \frac{n\sigmahat^2}{2} =
\frac{n\sigmahat^2}{2}\left(1 + \frac{2\theta}{n\sigmahat^2}\right),
\\ &
\tilde\mu_* = \frac{\mu_*+nv_*\ybar}{1+nv_*} =
\ybar\frac{1+\mu_*/(nv_*\ybar)}{1+1/(nv_*)}, \qquad
& &
\tilde v_* = \frac{v_*}{1+nv_*} =
\frac{1}{n}\frac{1}{1+1/(nv_*)}.
\end{alignedat}
$$


### μの周辺事前分布

共役事前分布における $\mu$ の周辺分布は次になる:

$$
\mu \sim
\mu_* + \sqrt{\frac{\theta}{\kappa}v_*}\;\op{TDist}(2\kappa).
$$

パラメータをBayes更新後に置き換えればこれは $\mu$ の周辺事後分布になる.


### 事前・事後予測分布

密度函数が

$$
p_*(y_\op{new}) =
\iint_{\R\times\R_{>0}}
p_N(y_\op{new}|\mu,v)p_{IGN}(\mu,v|\mu_*,v_*,\kappa,\theta)
\,d\mu\,dv
$$

で定義される事前予測分布は次になる:

$$
y_\op{new} \sim
\mu_* + \sqrt{\frac{\theta}{\kappa}(1+v_*)}\;\op{TDist}(2\kappa).
$$

パラメータをBayes更新後に置き換えればこれは事後予測分布になる.

<!-- #region -->
### Jeffreys事前分布の場合


正規分布モデルのJeffreys事前分布 $p_\op{Jeffreys}(\mu,v)$ は

$$
p_\op{Jeffreys}(\mu,v) \propto v^{-3/2}
$$

となることが知られている. ただし, これの $(\mu,v)\in\R\times\R_{>0}$ に関する積分は $\infty$ になるので, これはimproper事前分布である.

前節の逆ガンマ正規分布の密度函数と比較すると, Jeffreys事前分布に対応するパラメータ値は形式的に次になることがわかる:

$$
\kappa = 0, \quad
\theta = 0, \quad
v_* = \infty.
$$

そのとき, 前節のBayes更新の公式は次のようにシンプルになる:

$$
\tilde\kappa = \frac{n}{2}, \quad
\tilde\theta = \frac{n\sigmahat^2}{2}, \quad
\tilde\mu_* = \ybar, \quad
\tilde v_* = \frac{1}{n}.
$$

さらに, 前節の公式から, $n\to\infty$ のとき, 一般のパラメータ値に関するBayes更新の結果は, $n\to\infty$ のとき漸近的にこのJeffreys事前分布の場合に一致する.

さらに, Jeffreys事前分布の場合には

$$
\frac{\tilde\theta}{\tilde\kappa} = \sigmahat^2.
$$

ゆえに, $\mu$ に関する周辺事後分布は

$$
\mu \sim
\ybar + \frac{\sigmahat}{\sqrt{n}}\;\op{TDist}(n)
$$

になり, 事後予測分布は次になる:

$$
y_\op{new} \sim
\ybar + \sigmahat\sqrt{1+\frac{1}{n}}\;\op{TDist}(n).
$$
<!-- #endregion -->

### 通常の信頼区間と予測区間との比較

通常の $t$ 分布を使う平均の信頼区間と次の値の予測区間の構成では以下を使う:

$$
\frac{\ybar - \mu}{s\big/\!\sqrt{n}} \sim
\op{TDist}(n-1), \quad
\frac{y_\op{new} - \ybar}{s\sqrt{1+1/n}} \sim
\op{TDist}(n-1).
$$

ここで, $s^2$ は標本の不偏分散であり, $s$ はその平方根である:

$$
s^2 = \frac{1}{n-1}\sum_{i=1}^n(y_i - \ybar)^2 =
\frac{n\sigmahat^2}{n-1} > \sigmahat^2.
$$

したがって, 前節の結果と比較すると, Jeffreys事前分布の事後分布と予測分布による区間推定は, 以上の通常の区間推定よりも少し狭くなることがわかる.


### データの数値から事前分布を決めた場合

$a,b>0$ であると仮定する.

データの数値から共役事前分布のパラメータを次の条件によって決めたと仮定する:

$$
E[\mu] = \mu_* = \ybar, \quad
E[v] = \frac{\theta}{\kappa-1} = \sigmahat^2, \quad
\var(\mu) = v_* E[v] = a\sigmahat^2, \quad
\var(v) = \frac{E[v]^2}{\kappa-2} = b\sigmahat^4.
$$

これは次と同値である:

$$
\mu_* = \ybar, \quad
v_* = a, \quad
\kappa = 2 + \frac{1}{b}, \quad
\theta = \sigmahat^2\left(1 + \frac{1}{b}\right).
$$

これのBayes更新の結果は以下のようになる:

$$
\begin{alignedat}{2}
&
\tilde\kappa = 2 + \frac{1}{b} + \frac{n}{2} =
\frac{n}{2}\left(1 + \frac{2(2+1/b)}{n}\right)
& &
\to 2 + \frac{n}{2},
\\ &
\tilde\theta =
\sigmahat^2\left(1 + \frac{1}{b} + \frac{n}{2}\right) =
\frac{n\sigmahat^2}{2}\left(1 + \frac{2(1+1/b))}{n}\right)
& &
\to \sigmahat^2\left(1 + \frac{n}{2}\right),
\\ &
\tilde\mu_* = \frac{\ybar+nv_*\ybar}{1+nv_*} =
\ybar
& &
\to \ybar,
\\ &
\tilde v_* = \frac{a}{1+na} =
\frac{1}{n}\frac{1}{1+1/(na)}
& &
\to \frac{1}{n}.
\end{alignedat}
$$

以上における $\to$ は $a\to\infty$, $b\to\infty$ での極限を意味する.  その極限を取った後もJeffreys事前分布の場合とは一致しない.

以上の状況の下で,

$$
\frac{\tilde\theta}{\tilde\kappa} =
\sigmahat^2
\frac{1+2(1+1/b)/n}{1+2(2+1/b)/n} < \sigmahat^2,
\quad
\tilde v_* =
\frac{1}{n}\frac{1}{1+1/(na)} < \frac{1}{n}.
$$

ゆえに, この場合の区間推定はJeffreys事前分布の場合よりも狭くなる.

しかし, $n$ が大きければその違いは小さくなる.

```julia

```
